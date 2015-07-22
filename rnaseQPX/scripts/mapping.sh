#! /user/bin/bash

:'
this script accomplish 5 things:
1. map all paired end samples to reference woth bwa
2. sort the mapped contigs with samtools
3. remove duplicate contigs with picard
4. index contigs with samtools
5. count contigs with htseq
'

sample[1]=A1
sample[2]=A2
sample[3]=A3

ir=./trimmed/
dir=mapping5
ddir=rmdup5

extension=.trimmed.P.fastq.gz
reference=./genomeSRv015/QPX_v015.fasta
count=./genomeSRv015/QPX_v015.gff3

for i in 1 2 3
do
    sample=${sample[${i}]}
    bwa mem ${reference} \
        ${ir}${sample}R1${extension} \
        ${ir}${sample}R2${extension} | \
        samtools view -Shu - | \
        samtools sort - ./${dir}/${sample}.sorted

    samtools index ./${dir}/${sample}.sorted.bam

    htseq-count --format=bam \
        --stranded=no \
        --type=CDS --order=pos \
        --idattr=Name ./${dir}/${sample}.sorted.bam ${count} \
        > ./${ddir}/${sample}.htseq.counts.txt

    java -Xmx2g -jar ~/data/picard/picard.jar \
        MarkDuplicates \
        INPUT=./${dir}/${sample}.sorted.bam \
        OUTPUT=./${ddir}/${sample}.nodup.bam \
        METRICS_FILE=./${ddir}/${sample}.dup.metrics \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true

    rm -rf ./${dir}/${sample}.sorted.bam

    htseq-count --format=bam \
        --stranded=no \
        --type=CDS --order=pos \
        --idattr=Name ./${ddir}/${sample}.nodup.bam ${count} \
        > ./${ddir}/${sample}.htseq.counts.nodup.txt


done
