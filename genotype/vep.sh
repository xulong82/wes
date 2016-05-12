#!/bin/bash

module load vep

input="~/Dropbox/GitHub/wes/genotype/vep_input.txt"

cd /data/xwang/wes2

perl variant_effect_predictor.pl -i ${input} -o vep.txt
perl ${vep_path}/variant_effect_predictor.pl --offline -i ${input} -o vep.txt

module load perl
vep_path="/home/xwang/ensembl-tools-release-81/scripts/variant_effect_predictor"
cd /data/xwang/ADSP

par=$1

chr=chr$par
file=chr$par.vcf

# perl $path_vep/variant_effect_predictor.pl --offline --fork 20 \
#       --gmaf --symbol --protein --biotype --regulatory --force_overwrite \
#        -i ./VEP/$file -o ./VEP/$chr.txt

cd VEP
sed '/^##/d' $chr.txt > $chr.temp
mv $chr.temp $chr.txt

# awk '{$8 = $5"/"$6; $9 = "+"; print $2,$3,$3,$8,$9,$4}' ./IND0.01/$file > ./VEP/IND/$temp
#   perl $path_vep/variant_effect_predictor.pl --offline --fork 20 \
#        --gmaf --symbol --protein --biotype --regulatory --force_overwrite \
#        -i ./VEP/IND/$temp -o ./VEP/IND/$chr.txt

