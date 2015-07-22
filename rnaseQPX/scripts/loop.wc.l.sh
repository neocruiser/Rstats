#! /bin/bash

sample[1]=mmetsp0098
sample[2]=mmetsp001433
sample[3]=mmetsp00992
sample[4]=mmetsp001002
sample[5]=mmetsp0099
sample[6]=mmetsp00100

ir=/media/sf_docs/data/QPX-RNA-Seq/trimmed
dir=/media/sf_docs/data/mappingY4
ddir=/media/sf_docs/data/rmdupY4

counts=${ddir}/counts
realign=${ddir}/realign
call=${ddir}/call

for j in 1 2 3 4 5 6
do
    sample=${sample[${j}]}
    grep "N=2" ${call}/${sample}.last.call.vcf | wc -l

done
