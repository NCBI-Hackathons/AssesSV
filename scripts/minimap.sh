#!/bin/bash
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -q default
#$ -pe smp 16

minimap2 -t 16 \
        -a -k 19 \
        -O 5,56 -E 4,1 \
        -B 5 -z 400,50 -r 2k \
        --eqx --MD -Y \
        --secondary=no \
        -R "@RG\tID:${SAMPLE}_${FASTQ##*/}\tSM:${SAMPLE}" \
        "${REF}" "${FASTQ}" | \
        samtools sort -@ 16 - > "${OUTBAM}"
samtools index "${OUTBAM}"
