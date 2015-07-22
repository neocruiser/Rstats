#! /usr/bin/bash

zcat ./trim.A1R1.P.fastq.gz | fastqc ./trim.A1R1.P.fastq.gz && \
zcat ./trim.A1R2.P.fastq.gz | fastqc ./trim.A1R2.P.fastq.gz && \
zcat ./trim.A2R1.P.fastq.gz | fastqc ./trim.A2R1.P.fastq.gz && \
zcat ./trim.A2R2.P.fastq.gz | fastqc ./trim.A2R2.P.fastq.gz && \
zcat ./trim.A3R1.P.fastq.gz | fastqc ./trim.A3R1.P.fastq.gz && \
zcat ./trim.A3R2.P.fastq.gz | fastqc ./trim.A3R2.P.fastq.gz && \
zcat ./4.trim.A1R1.P.fastq.gz | fastqc ./4.trim.A1R1.P.fastq.gz && \
zcat ./4.trim.A1R2.P.fastq.gz | fastqc ./4.trim.A1R2.P.fastq.gz && \
zcat ./4.trim.A2R1.P.fastq.gz | fastqc ./4.trim.A2R1.P.fastq.gz && \
zcat ./4.trim.A2R2.P.fastq.gz | fastqc ./4.trim.A2R2.P.fastq.gz && \
zcat ./4.trim.A3R1.P.fastq.gz | fastqc ./4.trim.A3R1.P.fastq.gz && \
zcat ./4.trim.A3R2.P.fastq.gz | fastqc ./4.trim.A3R2.P.fastq.gz && \
zcat ./5.trim.A1R1.P.fastq.gz | fastqc ./5.trim.A1R1.P.fastq.gz && \
zcat ./5.trim.A1R2.P.fastq.gz | fastqc ./5.trim.A1R2.P.fastq.gz && \
zcat ./5.trim.A2R1.P.fastq.gz | fastqc ./5.trim.A2R1.P.fastq.gz && \
zcat ./5.trim.A2R2.P.fastq.gz | fastqc ./5.trim.A2R2.P.fastq.gz && \
zcat ./5.trim.A3R1.P.fastq.gz | fastqc ./5.trim.A3R1.P.fastq.gz && \
zcat ./5.trim.A3R2.P.fastq.gz | fastqc ./5.trim.A3R2.P.fastq.gz 

