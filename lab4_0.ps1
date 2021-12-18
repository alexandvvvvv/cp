##
# CPU: https://ark.intel.com/content/www/us/en/ark/products/201838/intel-core-i9-10980hk-processor-16m-cache-up-to-5-30-ghz.html
# VM with 4 cores
# VM with 4 GB RAM
# Ubuntu 21.04
# gcc version 10.3.0
##
$n0 = 2
$n1 = 700000
$n2 = 700000
$delta0 = ($n1 - $n0) / 100
$delta = ($n2 - $n1) / 10

$results_file = "lab4_results0_seq.txt"
$binary = "out/lab4_seq"
gcc -O3 -Wall -Werror main.c -lm -o $binary
Write-Output "Starting" | tee $results_file
for ($n = $n0; $n -le $n1; $n += $delta0) {
    $output = & $binary $n 4
    Write-Output $output[-1] | tee -a $results_file
}

$results_file = "lab4_results0_parall.txt"
$binary = "out/lab4"
gcc -O3 -Wall -Werror main.c -lm -fopenmp -o $binary
Write-Output "Starting" | tee $results_file
for ($n = $n0; $n -le $n1; $n += $delta0) {
    $output = & $binary $n 4
    Write-Output $output[-1] | tee -a $results_file
}
