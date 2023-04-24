module load plink/1.90b3

file=/path/snps.txt


cat $file | while read snp 
do
    /path/LD_eur_calculation.sh $snp
done
