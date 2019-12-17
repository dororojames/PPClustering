#include <stdio.h>
#include <stdlib.h>
#include "cluster.h"
#include <time.h>

int main(){
    srand(time(0));
    int N=10, dim=3;
    float data[N][dim];
    for (int i=0; i<N; i++) {
        for (int j=0; j<dim; j++) {
            data[i][j] = rand();
        }
    }
    clusters = Creat_Cluster()
}
