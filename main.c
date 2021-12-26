#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "./heapsort.c"
#include <float.h>
#include <unistd.h>
#include <pthread.h>
#include "./parallel.c"

#define A 432
#define ITERATIONS 1
unsigned int seed;

void log_array(char *message, double *array, int size)
{
  printf("\n%s", message);
  for (int i = 0; i < size; i++)
  {
    printf("%f ", array[i]);
  }
}

void fill_array(Chunk_t chunk)
{
  int N = chunk.offset + chunk.count;
  for (int i = chunk.offset; i < N; i++)
  {
    double value = (double)rand_r(&seed) / RAND_MAX;
    chunk.array1[i] = 1 + value * (A - 1);
  }
}

void fill_array2(Chunk_t chunk)
{
  int N = chunk.offset + chunk.count;
  for (int i = chunk.offset; i < N; i++)
  {
    double value = (double)rand_r(&seed) / RAND_MAX;
    chunk.array2[i] = A + value * (9 * A);
  }
}

void map_M1(double *array, int size)
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(array, size) schedule(guided, 4)
#endif
  for (int i = 0; i < size; i++)
  {
    array[i] = tanh(array[i]) - 1;
  }
}

void map_M2(double *array, int size)
{
  double *copy = calloc(size, sizeof(double));
  memcpy(copy + 1, array, sizeof(double) * (size - 1));
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

void merge(double *src_array, double *dest_array, int dest_size)
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(src_array, dest_array, dest_size) schedule(guided, 4)
#endif
  for (int i = 0; i < dest_size; i++)
  {
    dest_array[i] = fmax(src_array[i], dest_array[i]);
  }
}

double reduce(double *array, int size)
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
#pragma omp parallel for default(none) shared(array, size, min_value) reduction(+ \
                                                                                : result) schedule(guided, 4)
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

void sort(double *array, int size, int split_by_procs_num)
{
  int chunks = 2;

  //------------ calculate sub-arrays sizes ---------
  int *sizes = malloc(chunks * sizeof(int));

  int min_size = size / chunks;
  int inc_size = size % chunks;

  for (int i = 0; i < chunks; i++)
  {
    sizes[i] = min_size;

    if (inc_size > i)
      sizes[i] = sizes[i] + 1;
  }

  //------------- fill & sort sub-arrays -----------------
  double **temp = malloc(chunks * sizeof(double *));

  for (int i = 0; i < chunks; i++)
  {
    temp[i] = malloc(sizes[i] * sizeof(double));

    int start = 0;

    for (int j = 0; j < i; j++)
    {
      start += sizes[j];
    }

    for (int j = 0; j < sizes[i]; j++)
    {
      temp[i][j] = array[start + j];
    }

    heapSort(temp[i], sizes[i]);
  }

  //-------- use sub arrays to get sorted array --------
  int *visited_indecies = malloc(size * sizeof(int));

  for (int i = 0; i < size; i++)
  {
    int taken_from = 0;
    double next_value = DBL_MAX;

    for (int j = 0; j < chunks; j++)
    {
      int offset = visited_indecies[j];

      if (offset < sizeof(temp[j]))
      {
        if (next_value > temp[j][offset])
        {
          next_value = temp[j][offset];
          taken_from = j;
        }
      }
    }
    array[i] = next_value;
    visited_indecies[taken_from] = visited_indecies[taken_from] + 1;
  }

  free(sizes);
  for (int i = 0; i < chunks; i++)
  {
    free(temp[i]);
  }
  free(temp);
  free(visited_indecies);
}

#ifdef _OPENMP
void start_timer(double *T1)
{
  *T1 = omp_get_wtime();
}
#else
void start_timer(struct timeval *T1)
{
  gettimeofday(T1, NULL);
}
#endif

#ifdef _OPENMP
void stop_timer(double *T1)
{
  // double T2 = omp_get_wtime();
  // long delta_ms = 1000*(T2) - (*T1) * 1000;
  // printf("elapsed: %ld\n", delta_ms);
}
#else
void stop_timer(struct timeval *T1)
{
  // struct timeval T2;
  // gettimeofday(&T2, NULL);
  // long delta_ms = 1000*(T2.tv_sec - T1->tv_sec) + (T2.tv_usec - T1->tv_usec) / 1000;
  // printf("elapsed: %ld\n", delta_ms);
}
#endif

void *print_progress(void *i)
{
  while ((*(int *)i) < ITERATIONS)
  {
    double progress = (double)(*(int *)i) / ITERATIONS * 100;
    printf("Progress: %.2f%%\n", progress);
    sleep(1);
  }

  return 0;
}

int main(int argc, char *argv[])
{
  int N;
  unsigned int i;

  struct timeval T0, T1, T2;

  long delta_ms;
  N = atoi(argv[1]);

  gettimeofday(&T0, NULL);
  gettimeofday(&T1, NULL);
  int m_size = N;
  int m2_size = N / 2;

  double *M = malloc(m_size * sizeof(double));
  double *M2 = malloc(m2_size * sizeof(double));

  char *end;
  int num_threads = 1, is_parallel = 0, chunk_size = 4;
  if (argc > 2)
  {
    num_threads = (int)strtol(argv[2], &end, 10);
    if (num_threads > 1)
    {
      is_parallel = 1;
    }
    else
    {
      if (argc > 3)
      {
        chunk_size = (int)strtol(argv[3], &end, 10);
      }
    }
  }

  pthread_t progress_thread;
  pthread_create(&progress_thread, NULL, print_progress, &i);

  for (i = 0; i < ITERATIONS; i++)
  {
    seed = i;
    srand(seed * seed);

    //----------- Generate --------------//
    start_timer(&T1);
    if (is_parallel)
    {
      multiThreadComputing(M, M2, m_size, num_threads, fill_array, chunk_size);
      multiThreadComputing(M, M2, m2_size, num_threads, fill_array2, chunk_size);
    }
    else
    {
      fill_array(fullChunk(M, NULL, m_size));
      fill_array2(fullChunk(NULL, M2, m2_size));
    }
    stop_timer(&T1);

    //-------------- Map ----------------//
    start_timer(&T1);
    map_M1(M, m_size);
    map_M2(M2, m2_size);
    stop_timer(&T1);
    //log_array("Map M1: ", M, m_size);
    //log_array("Map M2: ", M2, m2_size);

    //------------- Merge ---------------//
    start_timer(&T1);
    merge(M, M2, m2_size);
    stop_timer(&T1);
    //log_array("Merge: ", M2, m2_size);

    //------------- Sort ----------------//
    start_timer(&T1);
    sort(M2, m2_size, 1);
    // heapSort(M2, m2_size);
    stop_timer(&T1);
    //log_array("Sort: ", M2, m2_size);

    //------------ Reduce ---------------//
    start_timer(&T1);
    reduce(M2, m2_size);
    stop_timer(&T1);
    // result = reduce(M2, m2_size);
  }
  free(M);
  free(M2);
  // printf("\nResult=%f", result);

  gettimeofday(&T2, NULL);
  pthread_join(progress_thread, NULL);
  delta_ms = 1000 * (T2.tv_sec - T0.tv_sec) + (T2.tv_usec - T0.tv_usec) / 1000;

  printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms);

  return 0;
}
