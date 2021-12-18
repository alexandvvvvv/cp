##
# CPU: https://ark.intel.com/content/www/us/en/ark/products/201838/intel-core-i9-10980hk-processor-16m-cache-up-to-5-30-ghz.html
# VM with 4 cores
# VM with 4 GB RAM
# Ubuntu 21.04
# gcc version 10.3.0
##
$nx = 266000
$n1 = 1900
$n2 = 700000
$minN = [math]::Min($nx / 2, $n1)
$delta = ($n2 - $minN) / 10

$binary = "out/lab4"
gcc -O3 -Wall -Werror main.c -lm -fopenmp -o $binary
for ($i = 1; $i -le 10; ++$i) {
    $results_file = "lab4_results1_$i.txt"
    Write-Output "Starting" | tee $results_file
    for ($n = $minN; $n -le $n2; $n += $delta) {
        $output = & $binary $n 4
        Write-Output $output[-1] | tee -a $results_file
    }    
}

