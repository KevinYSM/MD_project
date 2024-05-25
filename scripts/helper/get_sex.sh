 #!/bin/bash

 x_map=$(samtools idxstats $1 | grep "chrX\s" | cut -f 3)
 x_len=$(samtools idxstats $1 | grep "chrX\s" | cut -f 2)
 x_cov=$(echo "scale=5; ${x_map}/${x_len}" | bc)
echo "${x_map}.${x_len},${x_cov}"

 y_map=$(samtools idxstats $1 | grep "chrY\s" | cut -f 3)
 y_len=$(samtools idxstats $1 | grep "chrY\s" | cut -f 2)
 y_cov=$(echo "scale=5; ${y_map}/${y_len}" | bc)
 echo "${y_map}.${y_len}.${y_cov}"
 ratio=$(echo "scale=2; ${x_cov}/${y_cov}" | bc)
 echo $ratio

 if (( $(echo "$ratio > 4.00" | bc -l) )); then
     echo "F"
 else
     echo "M"
 fi