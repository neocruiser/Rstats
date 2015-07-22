#! /usr/bin/bash

dir=/media/sf_data/nodule/rmdup/
ddir=/home/neo/data/nodule/trinity

x=A1
y=A2
z=A3
b=A


    java -Xmx10g -jar /home/neo/data/picard/picard.jar \
        MergeSamFiles \
        I=${dir}${x}.nodup.bam \
        I=${dir}${y}.nodup.bam \
        I=${dir}${z}.nodup.bam \
        O=${dir}/${b}.bam \
        SO=coordinate \
        AS=true



/home/neo/data/QPX/trinityrnaseq/Trinity \
--genome_guided_bam ${dir}${b}.bam \
--genome_guided_max_intron 1000 \
--max_memory 10G \
--output ${ddir} \
--CPU 5

/home/neo/data/QPX/trinityrnaseq/Trinity --genome_guided_bam /media/sf_data/nodule/rmdup/A.bam --genome_guided_max_intron 1000 --max_memory 10G --output /home/neo/data/nodule/trinityA/ --CPU 8

