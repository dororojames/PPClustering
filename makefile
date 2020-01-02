all: *.c *.cu serial pthread omp cuda

serial: cluster.c
	gcc cluster.c -o cluster-s
pthread: cluster-pthread.c
	gcc -pthread -std=gnu99 -O2 -s cluster-pthread.c -o cluster-pthread
omp: cluster-omp.c
	gcc cluster-omp.c -o cluster-omp -fopenmp
cuda: cluster.cu
	nvcc cluster.cu -o cluster-cuda -L/usr/local/cuda/lib64 -lcurand
clean:
	rm cluster-s cluster-omp cluster-pthread cluster-cuda
