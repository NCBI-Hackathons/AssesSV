```bash
mkdir -p {subreads,10kb,15kb}/{fastx,minimap2,ngmlr}

## get fastas/fastqs
cd subreads/fastx
lftp -e "mget /giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/PacBio_fasta/*.fasta.gz" ftp://ftp-trace.ncbi.nlm.nih.gov
cd ../../

cd 10kb/fastx
lftp -e "get /giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_10kb/HG002.Q20.tar.bz2" ftp://ftp-trace.ncbi.nlm.nih.gov
tar jxvf PacBio_CCS_10kb/HG002.Q20.tar.bz2
cd ../../

cd 15kb/fastx
lftp -e "get /giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/HG002.Q20.tar.gz" ftp://ftp-trace.ncbi.nlm.nih.gov
tar gxvf PacBio_CCS_15kb/HG002.Q20.tar.bz2
cd ../../

## map fastas/fastqs with minimap2
cd subreads
parallel 'SAMPLE=HG002 REF=~/data/references/hs37d5.fa FASTQ={} OUTBAM={= s/subreads\.fasta\.gz/minimap2.bam/ =} qsub ~/scripts/minimap.sh' ::: fastx/*.fasta.gz
### wait
mv *.minimap2.bam minimap2
cd minimap2
ls *.bam > bam.list
samtools merge -@ 7 -b bam.list HG002.subreads.minimap2.hs37d5.bam
samtools index -@ 7 HG002.subreads.minimap2.hs37d5.bam
cd ../../

cd 10kb
parallel 'SAMPLE=HG002 REF=~/data/references/hs37d5.fa FASTQ={} OUTBAM={= s/Q20\.fastq/minimap2.bam/ =} qsub ~/scripts/minimap.sh' ::: fastx/*.fastq
### wait
mv *.minimap2.bam minimap2
cd minimap2
ls *.bam > bam.list
samtools merge -@ 7 -b bam.list HG002.10kb.minimap2.hs37d5.bam
samtools index -@ 7 HG002.10kb.minimap2.hs37d5.bam
cd ../../

cd 15kb
parallel 'SAMPLE=HG002 REF=~/data/references/hs37d5.fa FASTQ={} OUTBAM={= s/Q20\.fastq/minimap2.bam/ =} qsub ~/scripts/minimap.sh' ::: fastx/*.fastq
### wait
mv *.minimap2.bam minimap2
cd minimap2
ls *.bam > bam.list
samtools merge -@ 7 -b bam.list HG002.15kb.minimap2.hs37d5.bam
samtools index -@ 7 HG002.15kb.minimap2.hs37d5.bam
cd ../../

## map fastas/fastqs with ngmlr
cd subreads
parallel 'SAMPLE=HG002 REF=~/data/references/hs37d5.fa FASTQ={} OUTBAM={= s/subreads\.fasta\.gz/ngmlr.bam/ =} qsub ~/scripts/ngmlr.sh' ::: fastx/*.fasta.gz
### wait
mv *.ngmlr.bam ngmlr
cd ngmlr
ls *.bam > bam.list
samtools merge -@ 7 -b bam.list HG002.subreads.ngmlr.hs37d5.bam
samtools index -@ 7 HG002.subreads.ngmlr.hs37d5.bam
cd ../../

cd 10kb
parallel 'SAMPLE=HG002 REF=~/data/references/hs37d5.fa FASTQ={} OUTBAM={= s/Q20\.fastq/ngmlr.bam/ =} qsub ~/scripts/ngmlr.sh' ::: fastx/*.fastq
### wait
mv *.ngmlr.bam ngmlr
cd ngmlr
ls *.bam > bam.list
samtools merge -@ 7 -b bam.list HG002.10kb.ngmlr.hs37d5.bam
samtools index -@ 7 HG002.10kb.ngmlr.hs37d5.bam
cd ../../

cd 15kb
parallel 'SAMPLE=HG002 REF=~/data/references/hs37d5.fa FASTQ={} OUTBAM={= s/Q20\.fastq/ngmlr.bam/ =} qsub ~/scripts/ngmlr.sh' ::: fastx/*.fastq
### wait
mv *.ngmlr.bam ngmlr
cd ngmlr
ls *.bam > bam.list
samtools merge -@ 7 -b bam.list HG002.15kb.ngmlr.hs37d5.bam
samtools index -@ 7 HG002.15kb.ngmlr.hs37d5.bam
cd ../../
```
