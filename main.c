#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "./heapsort.c"

#define A 432

void log_array(char* message, double * array, int size)
{
  printf("\n%s", message);
  for (int i = 0; i < size; i++)
  {
    printf("%f ", array[i]);
  }
}

void fill_array(double * array, int size, unsigned int * seed, int min_value, int max_value) 
{
  for (int i = 0; i < size; i++) 
  {
    double value = (double) rand_r(seed) / RAND_MAX; 
    array[i] = min_value + value * (max_value - min_value);
  }
}

void map_M1(double * array, int size)
{
  for (int i = 0; i < size; i++)
  {
    array[i] = tanh(array[i]) - 1;
  }
}

void map_M2(double * array, int size)
{
  double * copy = malloc(size * sizeof(double));
  for (int i = 0; i < size; i++)
  {
    copy[i] = array[i];
  }
  
  array[0] = abs(cos(copy[0]));
  
  for (int i = 1; i < size; i++)
  {
    array[i] = fabs(cos(copy[i] + copy[i - 1]));
  }
  free(copy);
}

void merge(double * src_array, double * dest_array, int dest_size)
{
  for (int i = 0; i < dest_size; i++)
  {
    dest_array[i] = fmax(src_array[i], dest_array[i]);
  }
}

double reduce(double * array, int size)
{
  double min_value = 1.0;
  for (int i = 0; i < size; i++)
  {
    if (array[i] > 0)
    {
      min_value = array[i];
      break;
    }
  }
  
  double result = 0.0;
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
  int N;
  unsigned int i;
  struct timeval T1, T2;
  long delta_ms;
  N = atoi(argv[1]);
  gettimeofday(&T1, NULL); 
  int m_size = N;
  int m2_size = N / 2;
  double * M = malloc(m_size * sizeof(double));
  double * M2 = malloc(m2_size * sizeof(double));
  
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
    double result = reduce(M2, m2_size);
    printf("\nResult=%f", result);
  }
  
  gettimeofday(&T2, NULL); 
  delta_ms = 1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) / 1000;
  printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms);
  
  free(M);
  free(M2);
  
  return 0;
}
