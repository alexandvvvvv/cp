#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "./heapsort.c"
#include <float.h>
#include <unistd.h>

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

#ifdef _OPENMP
int omp_thread_count() {
    
    int n = 0;
    #pragma omp parallel reduction(+:n)
    n += 1;
    return n;
}
#endif


void fill_array(double * array, int size, int min_value, int max_value) 
{
  unsigned int seed = cos(size) * 10000;
  srand(seed);
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(array, size, min_value, max_value, seed) schedule(guided, 4)
  #endif
  for (int i = 0; i < size; i++) 
  {
    double value = (double) rand_r(&seed) / RAND_MAX; 
    array[i] = min_value + value * (max_value - min_value);
  }
}

void map_M1(double * array, int size)
{
#ifdef _OPENMP
  #pragma omp parallel for default(none) shared(array, size) schedule(guided, 4)
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
  #pragma omp parallel for default(none) shared(copy, array, size) schedule(guided, 4)
  #endif
  for (int i = 1; i < size - 1; i++)
  {
    copy[i] = array[i - 1];
  }
  
  array[0] = abs(cos(copy[0]));
  
#ifdef _OPENMP
  #pragma omp parallel for default(none) shared(copy, array, size) schedule(guided, 4)
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
  #pragma omp parallel for default(none) shared(src_array, dest_array, dest_size) schedule(guided, 4)
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
  #pragma omp parallel for default(none) shared(array, size, min_value) reduction(+: result) schedule(guided, 4)
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

void sort(double * array, int size, int split_by_procs_num) {
  int chunks = 2;
  
  #ifdef _OPENMP
  if (split_by_procs_num) {
    chunks = omp_get_num_procs();
  }
  #endif 
  
  //------------ calculate sub-arrays sizes ---------
  int * sizes = malloc(chunks * sizeof(int));
 
  int min_size = size / chunks;
  int inc_size = size % chunks;
  
  for (int i = 0; i < chunks; i++) {
     sizes[i] = min_size;
     
     if (inc_size > i) sizes[i] = sizes[i] + 1;
  }
  
  //------------- fill & sort sub-arrays -----------------
  double ** temp = malloc(chunks * sizeof(double *));

  #ifdef _OPENMP
  int initial_threads = omp_thread_count();
  omp_set_num_threads(chunks);
  #endif

#ifdef _OPENMP
  #pragma omp parallel for default(none) shared(chunks, sizes, temp, array) schedule(guided, 4)
  #endif
  for (int i = 0; i < chunks; i++) {
    temp[i] = malloc(sizes[i] * sizeof(double));
    
    int start = 0;
    
    for (int j = 0; j < i; j++) {
    	start += sizes[j];
    }
    
    for (int j = 0; j < sizes[i]; j++) {
       temp[i][j] = array[start + j];
    }
    
    heapSort(temp[i], sizes[i]);
    //log_array("sub array: ", temp[i], sizes[i]);
  }
  
  #ifdef _OPENMP
  omp_set_num_threads(initial_threads);
  #endif
  
  //-------- use sub arrays to get sorted array --------
  int * visited_indecies = malloc(size * sizeof(int));
  
  //do we need that ?
  for (int i = 0; i < chunks; i++) { visited_indecies[i] = 0; }
  
  for (int i = 0; i < size; i++) {
    int taken_from = 0;
    double next_value = DBL_MAX;
    
    for (int j = 0; j < chunks; j++) {
      if (visited_indecies[j] < sizes[j] && next_value > temp[j][visited_indecies[j]]) {
        next_value = temp[j][visited_indecies[j]];
        taken_from = j;
      }
    }
    //printf("next value = %f, taken_from = %d\n", next_value, taken_from);
    array[i] = next_value;
    visited_indecies[taken_from] = visited_indecies[taken_from] + 1;
  }  
}

int main(int argc, char* argv[])
{
  int N;
  unsigned int i;
  
  #ifdef _OPENMP
  double T1, T2;
  #else
  struct timeval T1, T2;
  #endif
  
  long delta_ms;
  N = atoi(argv[1]);
  
  #ifdef _OPENMP
  T1 = omp_get_wtime();
  #else
  gettimeofday(&T1, NULL); 
  #endif
  
  int m_size = N;
  int m2_size = N / 2;
  double * M = malloc(m_size * sizeof(double));
  double * M2 = malloc(m2_size * sizeof(double));
  #ifdef _OPENMP
  int threads_count = atoi(argv[2]);
  omp_set_num_threads(threads_count);
  #endif
  int iterations = 100;
#ifdef _OPENMP
  #pragma omp parallel
  #pragma omp sections
  {
  #pragma omp section
  {
  #endif
  for (i=0; i<iterations; i++) 
  {    
    //----------- Generate --------------//
    fill_array(M, m_size, 1, A);
    fill_array(M2, m2_size, A, 10 * A);

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
    sort(M2, m2_size, 0);
    
    //log_array("Sort: ", M2, m2_size);
    //------------ Reduce ---------------//
    reduce(M2, m2_size);
    
    //double result = reduce(M2, m2_size);
    //printf("\nResult=%f", result);
  }
  #ifdef _OPENMP
  T2 = omp_get_wtime();
  delta_ms = 1000*(T2) - (T1) * 1000;
  }
  #endif

  #ifdef _OPENMP
  #pragma omp section
  {
    while (i < iterations) {
      double progress = (double)i / iterations * 100;
      printf("Progress: %.2f\%%\n", progress);
      sleep(1);
    }
  }
  }
  #endif
  
  #ifndef _OPENMP
  gettimeofday(&T2, NULL); 
  delta_ms = 1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) / 1000;
  #endif
  
  printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms);
  
  free(M);
  free(M2);
  
  return 0;
}
