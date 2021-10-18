#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "./heapsort.c"
#include "fwBase.h"
#include "fwImage.h"
#include "fwSignal.h"
#define A 432

void log_array(char* message, float * array, int size)
{
  printf("\n%s", message);
  for (int i = 0; i < size; i++)
  {
    printf("%f ", array[i]);
  }
}

void fill_array(float * array, int size, unsigned int * seed, int min_value, int max_value) 
{
  for (int i = 0; i < size; i++) 
  {
    float value = (float) rand_r(seed) / RAND_MAX; 
    array[i] = min_value + value * (max_value - min_value);
  }
}

void map_M1(float * array, int size)
{
  fwsTanh_32f_A24(array, array, size);
  fwsAddC_32f(array, -1, array, size);
}

void map_M2(float * array, int size)
{
  float * copy = malloc(size * sizeof(float));
  copy[0] = 0;
  fwsCopy_32f(array, &copy[1], size - 1);

  fwsAdd_32f(array, copy, array, size);
  free(copy);

  fwsCos_32f_A24(array, array, size);
  fwsAbs_32f_I(array, size);
}

void merge(float * src_array, float * dest_array, int dest_size)
{
  fwsMaxEvery_32f_I(src_array, dest_array, dest_size);
}

float reduce(float * array, int size)
{
  float min_value = 1.0;
  for (int i = 0; i < size; i++)
  {
    if (array[i] > 0)
    {
      min_value = array[i];
      break;
    }
  }
  
  float result = 0.0;
  for (int i = 0; i < size; i++)
  {
    if ((int)(array[i] / min_value) % 2 == 0)
    {
       result += sin(array[i]);
    }
  }
  return result;
}

int main(int argc, char* argv[])
{
  unsigned int N;
  unsigned int i;
  struct timeval T1, T2;
  long delta_ms;
  N = atoi(argv[1]);
  gettimeofday(&T1, NULL); 
  int m_size = N;
  int m2_size = N / 2;
  float * M = malloc(m_size * sizeof(float));
  float * M2 = malloc(m2_size * sizeof(float));
  fwSetNumThreads(N);
  for (i=0; i<100; i++) 
  {
    unsigned int seed = i;
    srand(seed);
    
    //----------- Generate --------------//
    fill_array(M, m_size, &seed, 1, A);
    fill_array(M2, m2_size, &seed, A, 10 * A);

    //log_array("Initial M1: ", M, m_size);
    //log_array("Initial M2: ", M2, m2_size);
    //-------------- Map ----------------//
    map_M1(M, m_size);
    map_M2(M2, m2_size);
    
    //log_array("Map M1: ", M, m_size);
    //log_array("Map M2: ", M2, m2_size);
    //------------- Merge ---------------//
    merge(M, M2, m2_size);
    
    //log_array("Merge: ", M2, m2_size);
    //------------- Sort ----------------//
    heapSort(M2, m2_size);
    
    //log_array("Sort: ", M2, m2_size);
    //------------ Reduce ---------------//
    float result = reduce(M2, m2_size);
    printf("\nResult=%f", result);
  }
  
  gettimeofday(&T2, NULL); 
  delta_ms = 1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) / 1000;
  printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms);
  
  free(M);
  free(M2);
  
  return 0;
}
