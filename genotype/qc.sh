#!/bin/bash
#PBS -l mem=64gb,nodes=1:ppn=1,walltime=12:00:00

# Quality Flag 

cd /data/xwang/wes/qc

# files=`find /sdata/ADSP/processed_exome/casecontrol/*_SRR*/*_variants_filtered_highestsnpEff.vcf.gz`

files=`find /sdata/ADSP/processed_exome/enriched/*_SRR*/*_variants_filtered_highestsnpEff.vcf.gz`

for name1 in $files; do
    name2=`basename $name1`
    name3=${name2/_variants_filtered_highestsnpEff.vcf.gz/}

    echo $name3
    cp $name1 ./temp6.vcf.gz

    gzip -d temp6.vcf.gz # unzip the gz

    sed '/^#/d' temp6.vcf > temp1 # delete meta lines
    grep -v "LowCoverage" temp1 > temp2 # delete LowCoverage
    grep -v "VeryLowQual" temp2 > temp3 # delete VeryLowQual

    # take the snp id and genotype value columns
    awk '{print $1"\t"$2"\t"$7}' temp3 > temp4
    # header
    awk 'BEGIN {print "CHR POS QUAL"} {print $0}' temp4 > ${name3}.vcf

    rm temp[12346]*
done

