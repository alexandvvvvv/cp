#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "./heapsort.c"

#ifdef _OPENMP
#include <omp.h>
#endif

#define A 432

void log_array(char* message, double * array, int size)
{
  printf("\n%s", message);
  for (int i = 0; i < size; i++)
  {
    printf("%f ", array[i]);
  }
}

void fill_array(double * array, int size, unsigned int seed, int min_value, int max_value) 
{
  for (int i = 0; i < size; i++) 
  {
    double value = (double) rand_r(&seed) / RAND_MAX; 
    array[i] = min_value + value * (max_value - min_value);
  }
}

void map_M1(double * array, int size)
{
#ifdef _OPENMP
  #pragma omp parallel for default(none) shared(array, size)// schedule(runtime)
#endif
  for (int i = 0; i < size; i++)
  {
    array[i] = tanh(array[i]) - 1;
  }
}

void map_M2(double * array, int size)
{
  double * copy = malloc(size * sizeof(double));
  copy[0] = 0;
#ifdef _OPENMP
  #pragma omp parallel for default(none) shared(copy, array, size)// schedule(runtime)
#endif
  for (int i = 1; i < size - 1; i++)
  {
    copy[i] = array[i - 1];
  }
  
  array[0] = abs(cos(copy[0]));
  
#ifdef _OPENMP
  #pragma omp parallel for default(none) shared(copy, array, size)// schedule(runtime)
#endif
  for (int i = 1; i < size; i++)
  {
    array[i] = fabs(cos(array[i] + copy[i]));
  }
  free(copy);
}

void merge(double * src_array, double * dest_array, int dest_size)
{
#ifdef _OPENMP
  #pragma omp parallel for default(none) shared(src_array, dest_array, dest_size)// schedule(runtime)
#endif
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
#ifdef _OPENMP
  #pragma omp parallel for default(none) shared(array, size, min_value) reduction(+: result) //schedule(runtime)
#endif
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
  #ifdef _OPENMP
  int threads_count = atoi(argv[2]);
  omp_set_num_threads(threads_count);
  #if defined(SCHEDULE) && defined(CHUNKS)
  omp_set_schedule(SCHEDULE, CHUNKS);
  #endif
  #endif
  for (i=0; i<100; i++) 
  {
    unsigned int seed = i;
    srand(seed);
    
    //----------- Generate --------------//
    fill_array(M, m_size, seed, 1, A);
    fill_array(M2, m2_size, seed, A, 10 * A);

    // log_array("Initial M1: ", M, m_size);
    // log_array("Initial M2: ", M2, m2_size);
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
    reduce(M2, m2_size);
    // printf("\nResult=%f", result);
  }
  
  gettimeofday(&T2, NULL); 
  delta_ms = 1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) / 1000;
  printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms);
  
  free(M);
  free(M2);
  
  return 0;
}
