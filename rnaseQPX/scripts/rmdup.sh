#! /usr/bin/bash

count=./mmetsp0098/MMETSP0098.gff3

java -Xmx2g -jar ~/data/picard/picard.jar \
    MarkDuplicates \
        INPUT=./mapping3/A3.sorted.bam \
        OUTPUT=./rmdup3/A3.nodup.bam \
        METRICS_FILE=./rmdup3/A3.dup.metrics \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true

    #rm -rf ./${dir}/${sample}.sorted.bam

    htseq-count --format=bam \
        --stranded=no \
        --type=CDS --order=pos \
        --idattr=Name ./rmdup3/A3.nodup.bam ${count} \
        > ./rmdup3/A3.htseq.counts.nodup.txt
