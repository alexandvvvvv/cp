##
# CPU: https://ark.intel.com/content/www/us/en/ark/products/201838/intel-core-i9-10980hk-processor-16m-cache-up-to-5-30-ghz.html
# VM with 4 cores
# VM with 4 GB RAM
# Ubuntu 21.04
# gcc version 10.3.0
##
$n1 = 1900
$n2 = 700000
$delta = ($n2 - $n1) / 10

$threads_count = @(
    1,
    2,
    4,
    8,
    16
)

Write-Output "Compiling..."
gcc -O3 -Wall -Werror main.c -lm -o out/lab3_no_open_mp
gcc -O3 -Wall -Werror main.c -lm -fopenmp -o out/lab3_no_schedule

Write-Output "Checking compilation without OpenMP:"
$output = & out/lab3_no_open_mp $n1
Write-Output $output[-1]

Write-Output "N1: $n1, N2: $n2, D: $delta"

$results_file = 'lab3_results_no_schedule.txt'
Write-Output "Test without schedule" | tee $results_file

$binary = 'out/lab3_no_schedule'
foreach ($count in $threads_count) {
    Write-Output "Threads count: $count" | tee -a $results_file
    for ($n = $n1; $n -le $n2; $n += $delta) {
        $output = & $binary $n $count
        Write-Output $output[-1] | tee -a $results_file
    }
}

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
        $results_file = "lab3_results_${schedule}_${size}_chunks.txt"
        $binary = "out/lab3_schedule_${schedule}_${size}"
        gcc -O3 -Wall -Werror main.c -lm -fopenmp "-DSCHEDULE=omp_sched_$schedule" "-DCHUNKS=$size" -o $binary
        Write-Output "Schedule $schedule, chunk_size: $size" | tee $results_file
        for ($n = $n1; $n -le $n2; $n += $delta) {
            $output = & $binary $n $best_threads_count
            Write-Output $output[-1] | tee -a $results_file
        }
    }
}