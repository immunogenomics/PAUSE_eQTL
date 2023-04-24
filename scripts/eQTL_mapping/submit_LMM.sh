{
    local i=$1
        local j=$2

cat << EOF | bsub 

#!/bin/bash
#BSUB -J chr${i}_PC${j}
#BSUB -o /data/srlab2/qxiao/ITN/Joseph_scripts/Final/Qian_log/chr${i}_PC${j}.out
#BSUB -e /data/srlab2/qxiao/ITN/Joseph_scripts/Final/Qian_log/chr${i}_PC${j}.err
#BSUB -q big
#BSUB -M 50000
#BSUB -R 'rusage[mem=50000]'
#BSUB -n 4
#BSUB -N
#BSUB -u qxiao2@bwh.harvard.edu
#BSUB -R 'select[hname!=cn001]'
#BSUB -R 'select[hname!=cn002]'
#BSUB -R 'select[hname!=cn003]'
#BSUB -R 'select[hname!=cn004]'
#BSUB -R 'select[hname!=cn005]'

./LMM-logtpm.R ${i} ${j}

EOF
}

for i in {1..22}
#for i in {22..22}
do
 
  for j in 10 15 17 18 19 21 22 23 24 25 28 30 35
  
    do
      do_count $i $j
    done
done
