#!/bin/bash
# Author: XuLong Wang (xulong.wang@jax.org)

# Make TPED file

par=$1 # chromosome

vInf="meta_$par.txt"
geno="geno_$par.txt"

tped="chr$par.tped"

cd /data/xwang/wes/plink

nf=`awk 'NR==1 {print NF; exit}' ../golden/$geno`

awk '{print $1, $5, $6}' ../golden/$vInf > $par.0
awk '{$8 = 0} {print $2, $1, $8, $3}' ../golden/$vInf > $tped

for idx in `seq 1 $nf`; do
  awk -v col=$idx '{print $col}' ../golden/$geno > $par.1
  paste $par.0 $par.1 > $par.2
  awk '$4==0 {$5=$2; $6=$2} 
       $4==1 {$5=$2; $6=$3} 
       $4==2 {$5=$3; $6=$3} 
       {print $5, $6}' $par.2 > $par.3
  paste $tped $par.3 > $par.4
  mv $par.4 $tped
  rm $par.[123]
done

rm $par.0

