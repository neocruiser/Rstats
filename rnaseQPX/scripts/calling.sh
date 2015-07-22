#! /usr/bin/bash

:'
samtools -u for ouputing an uncompressed bcf file
-B : no baq computing for faster jobs
-d : depth of covreage, increase it to get precise depth of coverage
-f : decalre reference
-D : control the number of variant to keep per sample based on the depth of coverage
-C : reduce effect of reads with high mismatches
'

reference=./mmetsp0098/contigs.fa
sampleA1=./rmdup3/A1.nodup.bam
sampleA2=./rmdup3/A2.nodup.bam
sampleA3=./rmdup3/A3.nodup.bam

samtools mpileup -u -C50 -BQ0 -d10000000 -f ${reference} \
${sampleA1} ${sampleA2} ${sampleA3} |
bcftools view --min-ac 0 -g "^miss" - |
../bcftools-1.2/vcfutils.pl varFilter -D100 > ./calling3/mme98.var.vcf

# install vcf tools

