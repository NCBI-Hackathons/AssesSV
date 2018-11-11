/*
 * pipeline input parameters
 */
params.input_bam = "/home/ubuntu/data/HG002.10kb.Q20.GRCh38.pbmm2.bam"
params.outdir = "/home/ubuntu/pipeline-output"
params.reference_fasta = "/home/ubuntu/data/reference/"
params.bed = "/home/ubuntu/data/candidates.bam"

println "input bam: $params.input_bam"
println "outdir: $params.outdir"

ref = file(params.reference_fasta)
bed_path = file(params.smrtsv_bed)

bams = Channel.fromPath(params.input_bam)
bams.into { bam_sniffles; bam_pbsv; }

/*
 * sniffles
 */
process sniffles {
    publishDir "${params.outdir}/sniffles", mode:'copy'

    input:
    file bam from bam_sniffles

    output:
    file "sniffles.vcf" into sniffles_vcf_ch

    script:
    """

    """
}
/*
 * pbsv
 */
process pbsv {
    publishDir "${params.outdir}/pbsv", mode:'copy'

    input:
    file bam from bam_pbsv

    output:
    output:
    file "pbsv.svsig.gz" into pbsv_vcf_ch

    script:
    """
    pbsv discover $bam pbsv.svsig.gz
    """
}

/*
 * truevari
 */
process truevari {
    publishDir "${params.outdir}/truevari", mode:'copy'

    input:
    file pbsv_vcf from pbsv_vcf_ch
    file smartsv_vcf from smartsv_vcf_ch
    file sniffles_vcf from sniffles_vcf_ch

    output:
    file '*' into truevari_out_channel

    script:
    """

    """
}
