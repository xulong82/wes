#!/bin/bash

# Make compact VCF files

cd /data/xwang/wes/vcf

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
    awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$10}' temp3 > temp4
    # keep genotype info, discard others
    awk 'gsub(":.*", "", $6) {print $0}' temp4 > temp5
    # header
    awk 'BEGIN {print "CHR POS ID REF ALT GT"} {print $0}' temp5 > ${name3}.vcf

    rm temp[123456]*
done

