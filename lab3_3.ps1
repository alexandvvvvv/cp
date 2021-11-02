##
# CPU: https://ark.intel.com/content/www/us/en/ark/products/201838/intel-core-i9-10980hk-processor-16m-cache-up-to-5-30-ghz.html
# VM with 4 cores
# VM with 4 GB RAM
# Ubuntu 21.04
# gcc version 10.3.0
##
$best_schedule = 'guided'
$best_threads_count = 4
$best_chunk_size = 2
$n1 = 1900
$n2 = 700000
$delta = ($n2 - $n1) / 10

$optimizations = @(
    0,
    1,
    2,
    3
)

foreach ($level in $optimizations) {
    $results_file = "lab3_results_${best_schedule}_${best_chunk_size}_chunks_O${level}.txt"
    $binary = "out/lab3_schedule_${best_schedule}_${best_chunk_size}_O${level}"
    gcc "-O${level}" -Wall -Werror main.c -lm -fopenmp "-DSCHEDULE=omp_sched_$best_schedule" "-DCHUNKS=$best_chunk_size" -o $binary
    Write-Output "Schedule $best_schedule, chunk_size: $best_chunk_size, optimization: $level" | tee $results_file
    for ($n = $n1; $n -le $n2; $n += $delta) {
        $output = & $binary $n $best_threads_count
        Write-Output $output[-1] | tee -a $results_file
    }
}