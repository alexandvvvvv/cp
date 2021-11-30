##
# CPU: https://ark.intel.com/content/www/us/en/ark/products/201838/intel-core-i9-10980hk-processor-16m-cache-up-to-5-30-ghz.html
# VM with 4 cores
# VM with 4 GB RAM
# Ubuntu 21.04
# gcc version 10.3.0
##
$n1 = 1900
$best_threads_count = 4

$chunk_sizes = @(
    1,
    2,
    4,
    8
)
$schedules = @(
    'static',
    'dynamic',
    'guided'
)

foreach ($schedule in $schedules) {
    foreach ($size in $chunk_sizes) {
        $results_file = "lab3_results2_${schedule}_${size}_chunks.txt"
        $binary = "out/lab3_schedule_${schedule}_${size}"
        gcc -O3 -Wall -Werror main.c -lm -fopenmp "-DSCHEDULE=omp_sched_$schedule" "-DCHUNKS=$size" -o $binary
        Write-Output "Schedule $schedule, chunk_size: $size" | tee $results_file
        for ($n = 1; $n -le $n1; $n ++) {
            $output = & $binary $n $best_threads_count
            Write-Output $output[-1] | tee -a $results_file
        }
    }
}

$results_file = 'lab3_results2_no_open_mp.txt'
$binary = "out/lab3_no_open_mp"
gcc -O3 -Wall -Werror main.c -lm -o $binary
Write-Output "Schedule default, chunk_size: default" > $results_file
for ($n = 1; $n -le $n1; $n ++) {
    $output = & $binary $n $best_threads_count
    Write-Output $output[-1] >> $results_file
}