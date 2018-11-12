![alt text](https://github.com/NCBI-Hackathons/AssesSV/blob/master/graphics/assessv_v1.png)

## Assessing contribution of depth, read quality, and algorithm on Structural Variation calling.

Structural Variants (SVs), defined as insertions and deletions greater than 50 bp, inversions, and translocations, account for the majority of human genomic variation by total base pair count and have been implicated in genomic disorders. 

Pacific Biosciences has recently released two ~30-fold datasets of high-fidelity long reads, with mean insert sizes of 10 and 15 kb and read qualities above 99%. Because these high fidelity reads have a different error profile from PacBio subreads, it is unclear whether they will have similar performance in long read SV callers.  To test this, we cover the parameter space by comparing SV calls (using [truvari](https://github.com/spiralgenetics/truvari)) across error modes (high-fidelity vs raw subreads), read length (10 kbp vs 15 kbps vs 10-50 kbp ), aligner ([minimap2](https://github.com/lh3/minimap2) vs [NGM-LR](https://github.com/philres/ngmlr)), and variant caller ([Sniffles](https://github.com/fritzsedlazeck/Sniffles) vs [pbsv](https://github.com/PacificBiosciences/pbsv)).

## Datasets: 

[PacBio HG002 High-Fidelity QV20 10 kbp CCS, ~30-fold coverage](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_10kb/)
[PacBio HG002 High-Fidelity QV20 15 kbp CCS, ~30-fold coverage](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/)
[PacBio HG002 10-50 kbp Subreads, ~70-fold coverage](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/PacBio_fasta/)

[HG002 Tier 1 SV calls, version 0.6]()

## Conditions:

### Depth titration:
Starting at ~30-fold, downsample to 20-fold, 10-fold, 5-fold coverage

### Aligner:
#### NGMLR settings:
```bash
ngmlr --presets pacbio \
	--query "${FASTQ}" \
	--reference "${REF}" \
	--rg-id "${SAMPLE}_${FASTQ##*/}" \
	--rg-sm "${SAMPLE}"
```
#### minimap2 settings:
```bash
minimap2 -t 16 \
        -a -k 19 \
        -O 5,56 -E 4,1 \
        -B 5 -z 400,50 -r 2k \
        --eqx --MD -Y \
        --secondary=no \
        -R "@RG\tID:${SAMPLE}_${FASTQ##*/}\tSM:${SAMPLE}" \
        "${REF}" "${FASTQ}"
```
### Caller:
#### Sniffles settings:
```bash
sniffles -s 3 --skip_parameter_estimation -m !{bam} -v "!{basename}.!{depth}x.Sniffles_s3_ignoreParam.vcf"
```
#### pbsv settings:
```bash
pbsv discover !{bam} !{basename}.!{depth}x.pbsv.svsig.gz
pbsv call !{ref} !{basename}.!{depth}x.pbsv.svsig.gz !{basename}.!{depth}x.pbsv.vcf
```

### Workflow:
![alt text](https://github.com/NCBI-Hackathons/AssesSV/blob/master/graphics/Workflow_Full.png)
