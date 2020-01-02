#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <curand.h>

typedef struct group {
    // array of pointer of points
    int *p;

    // how many points in the group
    unsigned int length;
}Group;

typedef struct cluster {
    // array of pointer of groups
    Group** g;
    int num_point;
}Cluster;
int N;
int blockNum;

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

// CUDA
#define BLOCK_SIZE 1024

__constant__ int Num;

float *adj, *hadj;
int length;
int *gids, *rgids, *hgids;

int *gab, *hgab;


__device__ int index(const int i, const int j) {return (2*Num-i-1)*i/2+j-i-1;}
int hindex(const int i, const int j) {return (2*N-i-1)*i/2+j-i-1;}

void ShowAdj(const float* adj) {
    printf("============Adj============\n");
    for (int i = 0; i < N-1; i++) {
        for (int j = i+1; j < N; j++) {
            printf("%.4lf ", hadj[hindex(i, j)]);
        }
        printf("\n");
    }
    printf("===========================\n");
}

void ShowCluster(const Cluster* c, const int* gids) {
    printf("==========Cluster==========\n");
    printf("Groups=%d\n", length);
    for (int g = 0; g < length; g++) {
        int gid = hgids[g];
        printf("[");
        for (int i = 0; i < c->g[gid]->length; i++) {
            if (i) printf(" ");
            printf("%d", c->g[gid]->p[i]);
        }
        printf("]\n");
    }
    printf("===========================\n");
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

    cudaFree(gab);
    cudaFree(adj);
    cudaFree(gids);
    free(rgids);
    free(hgids);
    free(hgab);
    free(hadj);
}

void CheckError(const char name[]) {
    printf("%s\n", name);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        printf("Error: %s\n", cudaGetErrorString(err));
}

double randn() {
    int sign = (rand()%10)>5?1:-1;
    return rand()%10000000/(double)10000000.*sign;
}

double matmul(double input1[], double input2[], int dim) {
    double dis = 0;
    for (int i=0; i<dim; i++) {
        dis += input1[i] * input2[i];
    }
    return dis;
}

__global__ void Update(float *adj, int *gids, int *gab, int length) {
    int idx = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    if (idx >= length) return;
    int gid = gids[idx], ga = gab[0], gb = gab[1];
    if (gid < ga) {
        // printf("Update adj(%d, %d)\n", gid, ga);
        adj[index(gid, ga)] = max(adj[index(gid, ga)], adj[index(gid, gb)]);
    }
    else if (gid > ga && gid != gb) {
        // printf("Update adj(%d, %d)\n", gid, ga);
        int ra = min(gb, gid), rb = max(gb, gid);
        adj[index(ga, gid)] = max(adj[index(ga, gid)], adj[index(ra, rb)]);
    }
}

void Clusting() {
    double mind = 9999;
    for (int i = 0; i < length-1; i++) {
        for (int j = i+1; j < length; j++) {
            int gi = hgids[i], gj = hgids[j];
            double gd = hadj[hindex(gi, gj)];
            if (gd < mind) {
                mind = gd;
                hgab[0] = gi;
                hgab[1] = gj;
            }
        }
    }

    cudaMemcpy(gab, hgab, 2*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(gids, hgids, N*sizeof(int), cudaMemcpyHostToDevice);
    // printf("Merge %d %d\n", hgab[0], hgab[1]);


    int BN = length/BLOCK_SIZE;
    if (length % BLOCK_SIZE) BN += 1;
    dim3 blockNum(BN);
    Update<<<blockNum, BLOCK_SIZE>>>(adj, gids, gab, length);
    // CheckError("update");
    cudaDeviceSynchronize();
    // printf("Check U\n");
    cudaMemcpy(hadj, adj, N*(N-1)/2*sizeof(float), cudaMemcpyDeviceToHost);

    hgids[rgids[hgab[1]]] = hgids[length-1];
    rgids[hgids[length-1]] = rgids[hgab[1]];
    length -= 1;
}

int main(int argc, char *argv[]){
    curandGenerator_t gen;
    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(gen, 1234ULL);
    N = atoi(argv[1]);
    cudaMemcpyToSymbol(Num, &N, sizeof(int));

    /*
     adjacent list, which store the distance
     ex. adj[0][1] is the distance between point 0 and point 1
     */
    cudaMalloc((void**) &adj, N*(N-1)/2*sizeof(float));
    curandGenerateUniform(gen, adj, N*(N-1)/2);
    cudaMalloc((void**) &gids, N*sizeof(int));
    cudaMalloc((void**) &gab, 2*sizeof(int));
    hadj = (float*) malloc(N*(N-1)/2*sizeof(float));
    cudaMemcpy(hadj, adj, N*(N-1)/2*sizeof(float), cudaMemcpyDeviceToHost);
    rgids = (int*) malloc(N*sizeof(int));
    hgids = (int*) malloc(N*sizeof(int));
    hgab = (int*) malloc(2*sizeof(int));
    length = N;

    for (int i=0; i<N; ++i) {
      hgids[i] = i;
      rgids[i] = i;
    }

    Cluster *cluster = Create_Cluster(N);
    // ShowAdj(adj);
    // ShowCluster(cluster, gids);

    int cluster_size=1;
    while (length > cluster_size) {
        Clusting();

        Merge(cluster, hgab[0], hgab[1]);
        // ShowAdj(adj);
        // ShowCluster(cluster, gids);
    }

    // ShowCluster(cluster, gids);
    curandDestroyGenerator(gen);
    Clear(cluster);
}
