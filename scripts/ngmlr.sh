#!/bin/bash
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -q default
#$ -pe smp 16

ngmlr \
 --threads 16 \
 --presets pacbio \
 --reference "${REF}" \
 --rg-id "${SAMPLE}_${FASTQ##*/}"\
 --rg-sm "${SAMPLE}" \
 --query "${FASTQ}" | \
        samtools sort - > "${OUTBAM}"
samtools index "${OUTBAM}"
