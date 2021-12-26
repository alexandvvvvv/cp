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
double *copy;

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

void map_M1(Chunk_t chunk)
{
  for (int i = chunk.offset; i < chunk.offset + chunk.count; i++)
  {
    chunk.array1[i] = tanh(chunk.array1[i]) - 1;
  }
}

void map_M2(Chunk_t chunk)
{
  for (int i = chunk.offset; i < chunk.offset + chunk.count; i++)
  {
    chunk.array2[i] = fabs(cos(chunk.array2[i] + copy[i]));
  }
}

void merge(Chunk_t chunk)
{
  for (int i = chunk.offset; i < chunk.offset + chunk.count; i++)
  {
    chunk.array2[i] = fmax(chunk.array2[i], chunk.array1[i]);
  }
}

double reductionMin(double a, double b)
{
  if (isnan(b))
    return a;
  return MIN(a, b);
}

double reductionSum(double a, double b)
{
  if (isnan(b))
    return a;
  return a + b;
}

double get_min_positive(Chunk_t chunk)
{
  double result = DBL_MAX;
  double *M2 = chunk.array1;
  for (long j = chunk.offset; j < chunk.offset + chunk.count; j++)
  {
    if (M2[j] > 0 && M2[j] < result)
    {
      result = M2[j];
    }
  }
  return result;
}

double get_sinus_sum(Chunk_t chunk)
{
  double result = 0.0;
  double *M2 = chunk.array1;
  double min_value = chunk.array2[0];

  for (long j = chunk.offset; j < chunk.offset + chunk.count; j++)
  {
    if (((long)(M2[j] / min_value)) % 2 == 0)
    {
      result += sin(M2[j]);
    }
  }
  return result;
}

typedef struct
{
  double **temp;
  double *initial_array;
  double index;
  int *sizes;
} SortInfo_t;

void sort_part(void *args)
{
  SortInfo_t info = *(SortInfo_t *)args;
  int i = info.index;

  info.temp[i] = malloc(info.sizes[i] * sizeof(double));

  int start = 0;

  for (int j = 0; j < i; j++)
  {
    start += info.sizes[j];
  }

  for (int j = 0; j < info.sizes[i]; j++)
  {
    info.temp[i][j] = info.initial_array[start + j];
  }

  heapSort(info.temp[i], info.sizes[i]);
}

void sort(double *array, int size)
{
  int chunks = 4;

  //------------ calculate sub-arrays sizes ---------
  int *sizes = malloc(chunks * sizeof(int));

  int min_size = size / chunks;
  int inc_size = size % chunks;

  for (int i = 0; i < chunks; i++)
  {
    sizes[i] = min_size;

    if (inc_size > i)
    {
      sizes[i] = sizes[i] + 1;
    }
  }

  //------------- fill & sort sub-arrays -----------------
  double **temp = malloc(chunks * sizeof(double *));

  pthread_t *threads = malloc(chunks * sizeof(pthread_t));

  SortInfo_t *info = malloc(chunks * sizeof(SortInfo_t));
  for (int i = 0; i < chunks; i++)
  {
    info[i].initial_array = array;
    info[i].sizes = sizes;
    info[i].temp = temp;
    info[i].index = i;
    pthread_create(&threads[i], NULL, (void *)sort_part, (void *)&info[i]);
  }
  for (int i = 0; i < chunks; i++)
  {
    pthread_join(threads[i], NULL);
  }
  free(threads);
  free(info);

  //-------- use sub arrays to get sorted array --------
  int *visited_indecies = malloc(chunks * sizeof(int));

  for (int i = 0; i < chunks; i++)
  {
    visited_indecies[i] = 0;
  }

  for (int i = 0; i < size; i++)
  {
    int taken_from = 0;
    double next_value = DBL_MAX;

    for (int j = 0; j < chunks; j++)
    {
      int offset = visited_indecies[j];

      if (offset < sizes[j])
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

void start_timer(struct timeval *T1)
{
  gettimeofday(T1, NULL);
}
void stop_timer(struct timeval *T1)
{
  struct timeval T2;
  gettimeofday(&T2, NULL);
  long delta_ms = 1000*(T2.tv_sec - T1->tv_sec) + (T2.tv_usec - T1->tv_usec) / 1000;
  printf("elapsed: %ld\n", delta_ms);
}

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

  double result;
  for (i = 0; i < ITERATIONS; i++)
  {
    seed = i;
    srand(seed * seed);

    //----------- Generate --------------//
    start_timer(&T1);
    if (is_parallel)
    {
      compute_multi_thread(M, M2, m_size, num_threads, fill_array, chunk_size);
      compute_multi_thread(M, M2, m2_size, num_threads, fill_array2, chunk_size);
    }
    else
    {
      fill_array(get_full_chunk(M, NULL, m_size));
      fill_array2(get_full_chunk(NULL, M2, m2_size));
    }
    stop_timer(&T1);

    // log_array("Initial M1: ", M, m_size);
    // log_array("Initial M2: ", M2, m2_size);

    //-------------- Map ----------------//
    start_timer(&T1);
    copy = calloc(m2_size, sizeof(double));
    memcpy(copy + 1, M2, sizeof(double) * (m2_size - 1));
    copy[0] = 0;
    if (is_parallel)
    {
      compute_multi_thread(M, M2, m_size, num_threads, map_M1, chunk_size);
      compute_multi_thread(M, M2, m2_size, num_threads, map_M2, chunk_size);
    }
    else
    {
      map_M1(get_full_chunk(M, M2, m_size));
      map_M2(get_full_chunk(M, M2, m2_size));
    }
    free(copy);
    stop_timer(&T1);
    // log_array("Map M1: ", M, m_size);
    // log_array("Map M2: ", M2, m2_size);

    //------------- Merge ---------------//
    start_timer(&T1);
    if (is_parallel)
    {
      compute_multi_thread(M, M2, m2_size, num_threads, merge, chunk_size);
    }
    else
    {
      merge(get_full_chunk(M, M2, m2_size));
    }
    stop_timer(&T1);
    // log_array("Merge: ", M2, m2_size);

    //------------- Sort ----------------//
    start_timer(&T1);
    sort(M2, m2_size);
    // heapSort(M2, m2_size);
    stop_timer(&T1);
    // log_array("Sort: ", M2, m2_size);

    //------------ Reduce ---------------//
    start_timer(&T1);
    if (is_parallel)
    {
      double minPositive = reduce_multi_thread(M2, NULL, m2_size, num_threads, get_min_positive, reductionMin, DBL_MAX,
                                               chunk_size);
      result = reduce_multi_thread(M2, &minPositive, m2_size, num_threads, get_sinus_sum, reductionSum, 0,
                                   chunk_size);
    }
    else
    {
      double minPositive = get_min_positive(get_full_chunk(M2, NULL, m2_size));
      result = get_sinus_sum(get_full_chunk(M2, &minPositive, m2_size));
    }
    stop_timer(&T1);
  }
  free(M);
  free(M2);
  printf("\nResult=%f", result);

  gettimeofday(&T2, NULL);
  pthread_join(progress_thread, NULL);
  delta_ms = 1000 * (T2.tv_sec - T0.tv_sec) + (T2.tv_usec - T0.tv_usec) / 1000;

  printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms);

  return 0;
}
