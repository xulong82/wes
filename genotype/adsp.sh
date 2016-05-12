#!/bin/bash
#PBS -l mem=64gb,nodes=1:ppn=1,walltime=3:00:00

module load parallel
module load vcf-tools
module load bcftools
module load tabix

#--- Note: extract "BOTH_passed" variants: 30 mins

# dir="/sdata/ADSP/dbGaP-6143/files/vcf_snv"
# 
# cd /data/xwang/wes2/adsp
# files=`find ${dir}/c[12345]/*wes.chr${chr}.*.vcf.gz`
# 
# my_func1() {
#   cp $1 ./
#   name=`basename $1`
#   gzip -d ${name}
#   vcftools --vcf ${name/.gz}
#   vcftools --vcf ${name/.gz} --keep-filtered "BOTH_passed" --recode --out ${name/.vcf.gz}
#   bgzip ${name/.vcf.gz/.recode.vcf}
#   tabix -p vcf ${name/.vcf.gz/.recode.vcf.gz}
#   rm ${name/.gz}
# }
# 
# export -f my_func1
# 
# parallel my_func1 $1 ::: $files

#--- Note: merge c1-c5 per chromosome: 20 hours

# cd /data/xwang/wes2/adsp
# cs=`ls *.chr${chr}.*.vcf.gz`

# vcf-merge $cs | bgzip -c > ../merged/chr${chr}.merged.vcf.gz

#--- Note: efficient way to merge

# cd /data/xwang/wes2/adsp
# 
# file1=`ls c1_adsp.*chr${chr}.*.gz`
# output=chr${chr}.merged.vcf
# 
# cp $file1 ../merged2/${output}.gz
# cd ../merged2
# 
# gzip -d ${output}.gz
# grep -v "^##" ${output} > ${chr}_1
# mv ${chr}_1 ${output}
#   
# for idx in `seq 2 5`; do
#   file0=${file1/c1/c$idx}
#   cp ../adsp/${file0} ./
#   gzip -d ${file0}
#   grep -v "^##" ${file0/.gz} > ${chr}_1 && rm ${file0/.gz}
#   cut -d$'\t' -f10- ${chr}_1 > ${chr}_2 
#   paste ${output} ${chr}_2 > ${chr}_12
#   mv ${chr}_12 ${output} && rm ${chr}_[12]
# done

#--- Note: add file title

cd /data/xwang/wes2/merged2

sed '1d' chr${chr}.merged.filtered.recode.vcf > chr${chr}.tmp
cat ../header.txt chr${chr}.tmp > chr${chr}.merged.filtered.recode.vcf
rm chr${chr}.tmp

# sed -i '1s/^/##fileformat=VCFv4.2\n/' chr${chr}.merged.vcf
# sed '1s/4.2n/4.2\n/' chr${chr}.merged.vcf > ${chr}_1
# mv ${chr}_1 chr${chr}.merged.vcf

#--- Note: check sample order 

# cd /data/xwang/wes2/merged2
# ls *.vcf | parallel "bcftools query -l {} > {}.sp"
# ls *.sp | parallel diff chr1.merged.vcf.sp {}
# rm *.sp

#--- Note: filter variants

# cd /data/xwang/wes2/merged2

# name=chr${chr}.merged.vcf
# vcftools --vcf $name --maf 0.01 --recode --out ${name/.vcf/.filtered} 

#--- Note: concatenate vcf files

# cd /data/xwang/wes2/merged2

# files=`ls *.filtered.recode.vcf`
# vcf-concat $files > ../wes.vcf 

#--- Note: make plink files

# cd /data/xwang/wes2

# vcftools --vcf wes.vcf --plink-tped --out ./plink/wes 
 
