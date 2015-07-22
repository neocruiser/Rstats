#! /usr/bin/bash

~/data/QPX/trinityrnaseq/Trinity --seqType fq \
--left /media/Passport/MADL/QPX-RNA-Seq/MMETSP0098/MMETSP0098-Undescribed-sp--isolateNY0313808BC1.1.fastq.gz \
--right /media/Passport/MADL/QPX-RNA-Seq/MMETSP0098/MMETSP0098-Undescribed-sp--isolateNY0313808BC1.2.fastq.gz \
--quality_trimming_params "ILLUMINACLIP:~/data/Trimmomatic-0.33/adapters/TrueSeq3-PE-3.fa:2:30:10 TRAILING:3 MINLEN:36" \
--normalize_max_read_cov 50 \
--min_contig_length 200 \
--output ./assembl/trinity/ \
--max_memory 35G \
--CPU 10
