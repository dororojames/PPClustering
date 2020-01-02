#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>

typedef struct group Group;
typedef struct cluster Cluster;

struct group {
    // array of pointer of points
    int *p;

    // how many points in the group
    unsigned int length;
};

struct cluster {
    // array of pointer of groups
    Group** g;
    int num_point;
};

/*
    funct: create group and allocate the memory it needs
    parms:
        id: the id of group
    
    return:
        a pointer of Group
*/
Group* Create_Group(unsigned int id){
    Group* g = (Group*) malloc(sizeof(Group));
    g->p = (int*) malloc(1 * sizeof(int));
    g->length = 1;
    g->p[0] = id;

    return g;
}


/*
    funct: create cluster and allocate the memory it needs
    parms:
        n: the amount of points
    
    return:
        a pointer of Cluster
*/
Cluster* Create_Cluster(unsigned int n){
    Cluster* c = (Cluster*) malloc(sizeof(Cluster));
    c->g = (Group**) malloc(n * sizeof(Group*));

    c->num_point = n;

    int i;
    for(i = 0; i < n; ++i){
        c->g[i] = Create_Group(i);
    }

    return c;
}

/*
    funct: merge two groups
    parms:
        self: the group that absorb the other group
        other: the group that is going to be absorbed
    return:
        none
*/
void Merge_Group(Group* self, Group* other){
    unsigned int insert_position = self->length;
    self->length += other->length;
    self->p = (int*) realloc(self->p, self->length * sizeof(int));

    int i;
    for(i = 0; i < other->length; ++i){
        self->p[insert_position + i] = other->p[i];
    }

    free(other->p);
    free(other);
}

/*
    funct: merge two group in the cluster given by two group id
    parms:
        c: the cluster
        g1: group id of the group that will absorb the other group
        g2: group id of the group that will be absorbed
    return:
        none
*/
void Merge(Cluster* c, unsigned int g1, unsigned int g2){
    Merge_Group(c->g[g1], c->g[g2]);
    c->g[g2] = NULL;
}

void Clear(Cluster* c){
    int i;
    for(i = 0; i < c->num_point; ++i){
        if(c->g[i] == NULL){
            continue;
        }
        else{
            free(c->g[i]->p);
            free(c->g[i]);
        }
    }
}

unsigned int seed;
unsigned long long N, dim, length;
unsigned int num_threads;
unsigned int *tga, *tgb;
float *mind;
unsigned int ga, gb;
unsigned int *gids;
unsigned int *rgids;
    
int pthread_start[1000];
int pthread_num_point[1000];
int pthread_end[1000];

pthread_barrier_t t_barrier;

float** data;
float* adj;

float matmul(float input1[], float input2[], int dim) {
    float dis = 0;
    int i;
    
    for(i = 0; i < dim; ++i){
        dis += input1[i] * input2[i];
    }
    
    return dis;
}

int adjindex(const int i, const int j) {return (2*N-i-1)*i/2+j-i-1;}
float min(const float a, const float b) {return a<b?a:b;}
float max(const float a, const float b) {return a>b?a:b;}

void* initialize(void* threadId){
    int* tid = (int*)threadId;
    
    int i, j;
    
    // generate point position
    for(i = pthread_start[*tid]; i < pthread_end[*tid]; ++i){
        data[i] = (float*) malloc(dim * sizeof(float));
        
        for(j = 0; j < dim; ++j){
            int sign = (rand_r(&seed)%10)>5?1:-1;
            
            data[i][j] = rand_r(&seed)%10000000/(float)10000000.*sign;
        }
    }

    pthread_barrier_wait(&t_barrier);
    
    // initialize adj distance
    for(i = pthread_start[*tid]; i < pthread_end[*tid]; ++i){
        gids[i] = i;
        rgids[i] = i;
        
        for(j = i + 1; j < N; ++j){
            adj[adjindex(i, j)] = matmul(data[i], data[j], dim);
        }
    }
}


void* search_minimum(void* threadId){
    int* tid = (int*)threadId;
    
    int i, j;
    double local_min = 99999;
    
    for(i = pthread_start[*tid]; i < pthread_end[*tid]; ++i){
        for(j = i+1; j < length; ++j){
            int gi = gids[i], gj = gids[j];
            float gd = adj[adjindex(gi, gj)];
            
            if (gd < local_min) {
                
                local_min = gd;
                tga[*tid] = gi;
                tgb[*tid] = gj;
            }
        }
    }
    
    mind[*tid] = local_min;
}

void* update_adj(void* threadId){
    int* tid = (int*)threadId;
    
    int i;
    
    for(i = pthread_start[*tid]; i < pthread_end[*tid]; ++i){
        int gid = gids[i];
        
        if (gid < ga) {
            adj[adjindex(gid, ga)] = max(adj[adjindex(gid, ga)], adj[adjindex(gid, gb)]);
        }
        else if (gid > ga) {
            if (gid == gb){
                continue;
            }

            int ra = min(gb, gid), rb = max(gb, gid);
            adj[adjindex(ga, gid)] = max(adj[adjindex(ga, gid)], adj[adjindex(ra, rb)]);
        }
    }
}

void ShowAdj() {
    printf("============Adj============\n");
    for (int i = 0; i < N-1; i++) {
        for (int j = i+1; j < N; j++) {
            printf("%.4lf ", adj[adjindex(i, j)]);
        }
        printf("\n");
    }
    printf("===========================\n");
}

void ShowCluster(const Cluster* c) {
    printf("==========Cluster==========\n");
    printf("Groups=%lld\n", length);
    for (int g = 0; g < length; g++) {
        int gid = gids[g];
        printf("[");
        for (int i = 0; i < c->g[gid]->length; i++) {
            if (i) printf(" ");
            printf("%d", c->g[gid]->p[i]);
        }
        printf("]\n");
    }
    printf("===========================\n");
}

int main(int argc, char **argv)
{

    N = atoi(argv[1]);
    dim = atoi(argv[2]);
    
    unsigned int t;
    pthread_t* threads;
    
    num_threads = atoi(argv[3]);
    length = N;
    seed = time(NULL);
    
    threads = (pthread_t*) malloc(num_threads * sizeof(pthread_t));
    tga = (unsigned int*) malloc(num_threads * sizeof(unsigned int));
    tgb = (unsigned int*) malloc(num_threads * sizeof(unsigned int));
    mind = (float*) malloc(num_threads * sizeof(float));
    gids = (unsigned int*) malloc(N * sizeof(unsigned int));
    rgids = (unsigned int*) malloc(N * sizeof(unsigned int));
    adj = (float*) malloc(N * (N-1) / 2 * sizeof(float));
    
    data = (float **) malloc(N * sizeof(float*));

    int thread_id[num_threads];
    int rc;
    int point_per_thread;
    int working_thread;
    int remain_point;

    if(N < num_threads){
        point_per_thread = 1;
        working_thread = N;
        remain_point = 0;
    }
    else{
        point_per_thread = N / num_threads;
        working_thread = num_threads;
        remain_point = N % num_threads;
    }
    
    pthread_barrier_init (&t_barrier, NULL, working_thread);
    
    // initialize
    for (t = 0; t < working_thread; ++t) {
        thread_id[t] = t;
        
        if(t < remain_point){
            pthread_num_point[t] = point_per_thread + 1;
            pthread_start[t] = t * point_per_thread + t;
            
        }
        else{
            pthread_num_point[t] = point_per_thread;
            pthread_start[t] = t * point_per_thread + remain_point;
        }
        
        pthread_end[t] = pthread_start[t] + pthread_num_point[t];

        rc = pthread_create(&threads[t], NULL, initialize, (void*)&thread_id[t]);
        
        if (rc) {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }
    
    // initialize thread join
    for (t = 0; t < working_thread; ++t){
        pthread_join(threads[t], NULL);
    }
    
    Cluster *cluster = Create_Cluster(N);
    
    int cluster_size = 1;
    
    while (length > cluster_size) {
//        ShowAdj();
//        ShowCluster(cluster);
        
        float min_distance = 99999;

        // calculate range for every thread in searching minimun
        if(length < num_threads){
            point_per_thread = 1;
            working_thread = length;
            remain_point = 0;
        }
        else{
            point_per_thread = length / num_threads;
            working_thread = num_threads;
            remain_point = length % num_threads;
        }
        
        // search minimum
        for (t = 0; t < working_thread; ++t) {
            if(t < remain_point){
                pthread_num_point[t] = point_per_thread + 1;
                pthread_start[t] = t * point_per_thread + t;
                
            }
            else{
                pthread_num_point[t] = point_per_thread;
                pthread_start[t] = t * point_per_thread + remain_point;
            }
            
            pthread_end[t] = pthread_start[t] + pthread_num_point[t];
        
            rc = pthread_create(&threads[t], NULL, search_minimum, (void*)&thread_id[t]);
            
            if (rc) {
                printf("ERROR; return code from pthread_create() is %d\n", rc);
                exit(-1);
            }
        }
        
        // search minimum thread join
        for (t = 0; t < working_thread; ++t){
            pthread_join(threads[t], NULL);
            
            if(mind[t] < min_distance){
                min_distance = mind[t];
                ga = tga[t];
                gb = tgb[t];
            }
        }

        // update adj
        for (t = 0; t < working_thread; ++t) {
            rc = pthread_create(&threads[t], NULL, update_adj, (void*)&thread_id[t]);
            
            if (rc) {
                printf("ERROR; return code from pthread_create() is %d\n", rc);
                exit(-1);
            }
        }
        
        // update adj join
        for (t = 0; t < working_thread; ++t){
            pthread_join(threads[t], NULL);
        }

        gids[rgids[gb]] = gids[length-1];
        rgids[gids[length-1]] = rgids[gb];
        --length;
        
        Merge(cluster, ga, gb);
    }
    
//    ShowCluster(cluster);
    
    free(tga);
    free(tgb);
    free(gids);
    free(rgids);
    free(adj);
    
    for(int i = 0; i < N; ++i){
        free(data[i]);
    }
    free(data);
    
    pthread_barrier_destroy(&t_barrier);
    
    
    return 0;
}
