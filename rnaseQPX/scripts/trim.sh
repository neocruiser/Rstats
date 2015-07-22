#! /usr/bin/bash

:'
load program
load input files of reads
specify the names of output files (P: paired, U:Unpaired)
remove adapters
remove end of reads with low quality
remove reads less than 36
'
java -jar ~/data/Trimmomatic-0.33/trimmomatic-0.33.jar PE \
/media/Passport/MADL/Hard_clam_RNASeq/1_Index_1.A1_R1.fastq.gz \
/media/Passport/MADL/Hard_clam_RNASeq/HI.0615.001.Index_1.A1_R2.fastq.gz \
A1R1.trimmed.P.fastq.gz A1R1.trimmed.U.fastq.gz \
A1R2.trimmed.P.fastq.gz A1R2.trimmed.U.fastq.gz \
ILLUMINACLIP:TrueSeq3-PE-3.fa:2:30:10 \
TRAILING:3 \
MINLEN:36
