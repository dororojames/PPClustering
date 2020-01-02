#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

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
int N, Nthrds;

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


int adjindex(const int i, const int j) {return (2*N-i-1)*i/2+j-i-1;}

void ShowAdj(const float* adj, const int n) {
    printf("============Adj============\n");
    for (int i = 0; i < n-1; i++) {
        for (int j = i+1; j < n; j++) {
            printf("%.4lf ", adj[adjindex(i, j)]);
        }
        printf("\n");
    }
    printf("===========================\n");
}

void ShowCluster(const Cluster* c, const int* gids, const int length) {
    printf("==========Cluster==========\n");
    printf("Groups=%d\n", length);
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

void Clear(Cluster* c, int* gids, int* rgids, float* adj){
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
    free(gids);
    free(rgids);
    free(adj);
}

float randn() {
    int sign = (rand()%10)>5?1:-1;
    return rand()%10000000/(float)10000000.*sign;
}

float min(const float a, const float b) {return a<b?a:b;}
float max(const float a, const float b) {return a>b?a:b;}

float matmul(float input1[], float input2[], int dim) {
    float dis = 0;
//    #pragma omp parallel for reduction(+:dis)
    for (int i=0; i < dim; ++i) {
        dis += input1[i] * input2[i];
    }
    return dis;
}


void FindNearest(float* adj, int* gids, int* rgids, int* length, int *ga, int *gb) {
    float mind = 9999;
    int L = *length;
    #pragma omp parallel
    {
        int id, istart, iend, lga, lgb;
        float lmind = 9999;
        id = omp_get_thread_num();
        istart = id * L / Nthrds;
        iend = (id+1) * L / Nthrds;
        if (id == Nthrds-1) iend = L;
        for (int i = istart; i < iend; ++i) {
            for (int j = i+1; j < L; ++j) {
                int gi = gids[i], gj = gids[j];
                float gd = adj[adjindex(gi, gj)];
                if (gd < lmind) {
    //                printf("{%d %d}", gi, gj);
                    lmind = gd;
                    lga = gi;
                    lgb = gj;
                }
            }
        }
        
        # pragma omp barrier
        if (lmind < mind) {
            mind = lmind;
            *ga = lga;
            *gb = lgb;
        }

//        printf("Merge %d %d\n", *ga, *gb);
        
        for (int i = istart; i < iend; ++i) {
            int gid = gids[i];
            if (gid < *ga) {
    //            printf("Update adj(%d, %d)\n", gid, *ga);
                adj[adjindex(gid, *ga)] = max(adj[adjindex(gid, *ga)], adj[adjindex(gid, *gb)]);
            }
            else if (gid > *ga && gid != *gb) {
    //            printf("Update adj(%d, %d)\n", gid, *ga);
                int ra = min(*gb, gid), rb = max(*gb, gid);
                adj[adjindex(*ga, gid)] = max(adj[adjindex(*ga, gid)], adj[adjindex(ra, rb)]);
            }
        }
    }

    
//    #pragma omp parallel
//    {
//        int id, istart, iend;
//        id = omp_get_thread_num();
//        istart = id * length / Nthrds;
//        iend = (id+1) * length / Nthrds;
//        if (id == Nthrds-1) iend = length;
//
//    }

    
    gids[rgids[*gb]] = gids[L-1];
    rgids[gids[L-1]] = rgids[*gb];
    *length -= 1;
}

int main(int argc, char *argv[]){
    srand(66); //time(0)
    N = atoi(argv[1]);
    int dim = atoi(argv[2]);
    Nthrds = atoi(argv[3]);
    omp_set_num_threads(Nthrds);
    float data[N][dim];
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < dim; j++) {
            data[i][j] = randn();
        }
    }
    
//    for (int i = 0; i < N; i++) {
//        for (int j = 0; j < dim; j++) {
//            printf("%lf ", data[i][j]);
//        }
//        printf("\n");
//    }
    
    /*
        adjacent list, which store the distance
        ex. adj[0][1] is the distance between point 0 and point 1
    */
    float *adj = (float*) malloc(N*(N-1)/2 * sizeof(float));
    
    int length = N;
    int *gids = (int*) malloc(N*sizeof(int));
    int *rgids = (int*) malloc(N*sizeof(int));
    
    for (int i = 0; i < N-1; ++i) {
        for (int j = i+1; j < N; ++j) {
            float d = matmul(data[i], data[j], dim=dim);
            adj[adjindex(i, j)] = d;
        }
    }
    for (int i=0; i<N; ++i) {
        gids[i] = i;
        rgids[i] = i;
    }
    
    Cluster *cluster = Create_Cluster(N);
//    ShowAdj(adj, N);
//    ShowCluster(cluster, gids, length);
    
    int cluster_size=1;
    while (length > cluster_size) {
        int ga = 0, gb = 0;
        FindNearest(adj, gids, rgids, &length, &ga, &gb);
        
        Merge(cluster, ga, gb);
//        ShowAdj(adj, N);
//        ShowCluster(cluster, gids, length);
    }
//    ShowCluster(cluster, gids, length);
    Clear(cluster, gids, rgids, adj);
}
