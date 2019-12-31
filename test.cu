#include <cuda.h>
#include <stdio.h>
#include <time.h>


#define tbp 512
#define nblocks 1
__global__ void kernel_min(int *a, int *d)
{
  __shared__ int sdata[tbp]; //"static" shared memory

  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  sdata[tid] = a[i];

  __syncthreads();
  for(unsigned int s=tbp/2 ; s >= 1 ; s=s/2)
  {
    if(tid < s)
    {
      if(sdata[tid] > sdata[tid + s])
      {
        sdata[tid] = sdata[tid + s];
      }
    }
    __syncthreads();
  }
  if(tid == 0 )
  {
    d[blockIdx.x] = sdata[0];
  }
}

int main()
{
  int i;
  const int N=tbp*nblocks;
  srand(time(NULL));

  int *a;
  a = (int*)malloc(N * sizeof(int));
  int *d;
  d = (int*)malloc(nblocks * sizeof(int));

  int *dev_a, *dev_d;

  cudaMalloc((void **) &dev_a, N*sizeof(int));
  cudaMalloc((void **) &dev_d, nblocks*sizeof(int));
  int mmm=9999;
  for( i = 0 ; i < N ; i++)
  {
    a[i] = rand()%100/100;
    //printf("%d",a[i]);
    if(mmm>a[i]) mmm=a[i];

  }
  printf("");
  printf("");
  printf("");
  printf("");
  cudaMemcpy(dev_a , a, N*sizeof(int),cudaMemcpyHostToDevice);

  kernel_min <<<nblocks, tbp>>>(dev_a,dev_d);

  cudaMemcpy(d, dev_d, nblocks*sizeof(int),cudaMemcpyDeviceToHost);

  printf("cpu min %d, gpu_min = %d\n",mmm,d[0]);

  cudaFree(dev_a);
  cudaFree(dev_d);

  printf("");

  return 0;
}
