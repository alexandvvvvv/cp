##
# CPU: https://ark.intel.com/content/www/us/en/ark/products/201838/intel-core-i9-10980hk-processor-16m-cache-up-to-5-30-ghz.html
# VM with 4 cores
# VM with 4 GB RAM
# Ubuntu 21.04
# gcc version 10.3.0
##
$n1 = 1900
$n2 = 700000
$delta8 = ($n2 - $n1) / 10
$delta9 = ($n2 - $n1) / 4

Write-Output "N1: $n1, N2: $n2, D8: $delta8, D9: $delta9"

$binaries = @(
    './out/lab1-seq',
    './out/lab1-par-1',
    './out/lab1-par-2',
    './out/lab1-par-4',
    './out/lab1-par-8',
    './out/lab1-par-16'
)

Write-Output 'point 8'

foreach ($binary in $binaries) {
    Write-Output $binary
    for ($n = $n1; $n -le $n2; $n += $delta8) {
        $output = & $binary $n
        Write-Output $output[-1]   
    }
}

Write-Output ''
Write-Output 'point 9'

foreach ($binary in $binaries) {
    Write-Output $binary
    for ($n = $n1; $n -le $n2; $n += $delta9) {
        $output = & $binary $n
        Write-Output $output[-2]
        Write-Output $output[-1]   
    }
}