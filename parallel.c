#include <pthread.h>

typedef struct
{
    double *array1;
    double *array2;
    long offset;
    long count;
} Chunk_t;

typedef struct
{
    int threadId;
    void *func;
} ThreadInfo_t;

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

pthread_mutex_t mutex;
double *array1;
double *array2;
long arrayLength;
long currentOffset;
long chunkSize;

double reductionResult;
double (*reductionFunction)(double, double);

Chunk_t getNextChunk(int threadId, double intermediateReductionResult)
{
    pthread_mutex_lock(&mutex);
    if (reductionFunction != NULL)
    {
        reductionResult = (*reductionFunction)(reductionResult, intermediateReductionResult);
    }
    Chunk_t ans;
    ans.array1 = array1;
    ans.array2 = array2;
    
    if (currentOffset >= arrayLength)
    {
        ans.count = 0;
        ans.offset = 0;
    }
    else
    {
        ans.offset = currentOffset;
        ans.count = MIN(chunkSize, arrayLength - currentOffset);
    }
    currentOffset = currentOffset + chunkSize;
    pthread_mutex_unlock(&mutex);
    return ans;
}

Chunk_t fullChunk(double *_array1, double *_array2, long size)
{
    Chunk_t chunk;
    chunk.array1 = _array1;
    chunk.array2 = _array2;
    chunk.count = size;
    chunk.offset = 0;
    return chunk;
}

void threadFunc(void *args)
{
    ThreadInfo_t info = *(ThreadInfo_t *)args;
    Chunk_t chunk = getNextChunk(info.threadId, NAN);
    double (*f)(Chunk_t);
    f = info.func;
    while (chunk.count > 0)
    {
        double intermediateRes = f(chunk);
        chunk = getNextChunk(info.threadId, intermediateRes);
    }
}

double multiThreadComputingReduction(double *_array1, double *_array2, long _arrayLength, int threadsNum, void *func,
                                     void *_reductionFunction, double reductionInitValue, int _chunkSize)
{
    pthread_mutex_init(&mutex, NULL);

    pthread_t *thread = malloc(threadsNum * sizeof(pthread_t));
    ThreadInfo_t *info = malloc(threadsNum * sizeof(ThreadInfo_t));

    array1 = _array1;
    array2 = _array2;
    arrayLength = _arrayLength;
    chunkSize = _chunkSize;
    currentOffset = 0;
    reductionResult = reductionInitValue;
    reductionFunction = _reductionFunction;

    for (int j = 0; j < threadsNum; j++)
    {
        info[j].threadId = j + 1;
        info[j].func = func;
        pthread_create(&thread[j], NULL, (void *)threadFunc, (void *)(info + j));
    }

    for (int i = 0; i < threadsNum; i++)
    {
        pthread_join(thread[i], NULL);
    }

    free(thread);
    free(info);
    pthread_mutex_destroy(&mutex);

    return reductionResult;
}

double multiThreadComputing(double *_array1, double *_array2, long _arrayLength, int threadsNum, void *func,
                            int _chunkSize)
{
    return multiThreadComputingReduction(_array1, _array2, _arrayLength, threadsNum, func,
                                         NULL, NAN, _chunkSize);
}