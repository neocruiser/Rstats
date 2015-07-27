#! /user/bin/bash

:'
this script accomplish 6 things:
0. trimm reads
1. map all paired end samples to reference with bwa
2. sort the mapped reads with picard
3. remove duplicate reads with picard
4. index reads with samtools
5. count reads with htseq
'

sample[1]=C1
sample[2]=C2
sample[3]=C3

ir=/media/sf_data/nodule
irr=/media/sf_data/nodule/trimmed

dir=/media/sf_data/nodule/mapping
ddir=/media/sf_data/nodule/rmdup

extension=.trimmed.P.fastq.gz
reference=/media/sf_data/genomeSRv015/QPX_v015.fasta

## essential for calling SNPs
RG[1]='@RG\tID:noduleC1\tSM:MA\tPL:illumina\tLB:noduleC1\tPU:genomeSRv015'
RG[2]='@RG\tID:noduleC2\tSM:MA\tPL:illumina\tLB:noduleC2\tPU:genomeSRv015'
RG[3]='@RG\tID:noduleC3\tSM:MA\tPL:illumina\tLB:noduleC3\tPU:genomeSRv015'

stats=/media/sf_data/nodule/stats
count=/media/sf_data/genomeSRv015/QPX_v015.gff3


for i in 1 2 3
do
    sample=${sample[${i}]}
    java -Xmx10g -jar /home/neo/data/Trimmomatic-0.33/trimmomatic-0.33.jar PE \
        ${ir}/${sample}R1.fastq.gz \
        ${ir}/${sample}R2.fastq.gz \
        ${irr}/${sample}.1.trimmed.P.fastq.gz \
        ${irr}/${sample}.1.trimmed.U.fastq.gz \
        ${irr}/${sample}.2.trimmed.P.fastq.gz \
        ${irr}/${sample}.2.trimmed.U.fastq.gz \
        ILLUMINACLIP:TrueSeq3-PE-3.fa:2:30:10 \
        SLIDINGWINDOW:4:15 \
        TRAILING:5 \
        CROP:70 \
        MINLEN:30

    rm -f ${irr}/${sample}.1.trimmed.U.fastq.gz
    rm -f ${irr}/${sample}.2.trimmed.U.fastq.gz

done


## create dictionary and index of reference
    java -jar /home/neo/data/picard/picard.jar \
        CreateSequenceDictionary \
        R=${reference} \
        O=/media/sf_data/genomeSRv015/QPX_v015.dict

    samtools faidx ${reference}

## Map | Sort | remove duplicates
for i in 1 2 3
do
    sample=${sample[${i}]}
    RG=${RG[${i}]}
    bwa mem -M \
        -R ${RG} \
        -p ${reference} \
        ${irr}/${sample}.1${extension} \
        ${irr}/${sample}.2${extension} \
    > ${dir}/${sample}.sam

    java -Xmx10g -jar /home/neo/data/picard/picard.jar \
        SortSam \
        INPUT=${dir}/${sample}.sam \
        OUTPUT=${dir}/${sample}.sorted.bam \
        SORT_ORDER=coordinate

    java -Xmx10g -jar /home/neo/data/picard/picard.jar \
        MarkDuplicates \
        INPUT=${dir}/${sample}.sorted.bam \
        OUTPUT=${ddir}/${sample}.nodup.bam \
        METRICS_FILE=${stats}/${sample}.dup.metrics \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true

    htseq-count --format=bam \
        --stranded=no \
        --type=CDS --order=pos \
        --idattr=Name ${ddir}/${sample}.nodup.bam ${count} \
        > ${stats}/${sample}.htseq.counts.nodup.txt

done

