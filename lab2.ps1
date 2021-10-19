##
# CPU: https://ark.intel.com/content/www/us/en/ark/products/201838/intel-core-i9-10980hk-processor-16m-cache-up-to-5-30-ghz.html
# VM with 4 cores
# VM with 4 GB RAM
# Ubuntu 21.04
# gcc version 10.3.0
##
$n1 = 1900
$n2 = 700000
$k = 4
$delta8 = ($n2 - $n1) / 10
$delta9 = ($n2 - $n1) / 4

Write-Output "N1: $n1, N2: $n2, K: $k, D8: $delta8, D9: $delta9" | tee ./lab2_results

for ($m = 1; $m -le $k * 2; $m++) {
    Write-Output "M: $m" | tee -a ./lab2_results
    for ($n = $n1; $n -le $n2; $n += $delta8) {
        $output = ./out/lab2 $n $m
        Write-Output $output[-1] | tee -a ./lab2_results
    }
}

# Write-Output "Prepare to scan"
# Read-Host

# for ($m = 1; $m -le $k * 2; $m++) {
#     Write-Output "M: $m"
#     $output = ./out/lab2 $n2 $m
#     Write-Output $output[-1]
# }
