#!/bin/bash

function singleread(){
    for i in $1_*_R1_001.fastq
    do
    cat $i >> $1_R1.fastq
    done
    echo "Mapping $1; $2 cores"
    bowtie2 -x TAIR10 -U $1_R1.fastq -S $1.sam -p $2
    samtools view -@ $2 -bS $1.sam > $1.bam
    samtools sort -@ $2 $1.bam > $1.sort.bam
    samtools index -@ $2 $1.sort.bam
}

# 引数の取得
name1="$1"
cores="$2"

singleread "$name1" "$cores"