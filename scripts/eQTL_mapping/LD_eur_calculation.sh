#!/bin/sh
snp=$1
chr=$( echo $snp | cut -d "_" -f1 )
OFILE=/path/$snp.ld
rm -f $OFILE

echo $chr

N=$( cat /path/reference/1KG/EUR_plink/chr$chr.m2.M2.nodupid.bim |
   grep -F -w $snp - | wc -l )

if [ $N -gt 0 ] ; then
   plink \
   --bfile  /path/1KG/EUR_plink/chr$chr.m2.M2.nodupid  \
   --r2  \
   --allow-no-sex   \
   --threads 1  \
   --memory  2000   \
   --ld-window 500000   \
   --ld-window-r2 0.8 \
   --ld-snp  $snp  \
   --silent \
   --out  $OFILE.tmp01
   
   cat $OFILE.tmp01.ld |
   awk '{print $6,$7}' > $OFILE
   
else
   echo "SNP_B R2" > $OFILE
   echo "NA NA" >> $OFILE
fi

gzip -f $OFILE
rm -f $OFILE.tmp01*

