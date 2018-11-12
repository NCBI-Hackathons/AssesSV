nextflow run -resume  https://github.com/NCBI-Hackathons/AssesSV \
	-with-report -with-dag \
	--input_bams '*.bam' \
	--truth_vcf ~/data/giab-SV/chr12.HG002_SVs_Tier1_v0.6.vcf.gz \
	--bed ~/data/giab-SV/chr12.HG002_SVs_Tier1_v0.6.bed

