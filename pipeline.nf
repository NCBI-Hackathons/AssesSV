/*
 * pipeline input parameters
 */
params.input_bam = "/home/ubuntu/data/HG002.10kb.Q20.GRCh38.pbmm2.bam"
params.outdir = "/home/ubuntu/pipeline-output"
params.reference_fasta = "/home/ubuntu/data/reference/"
params.bed = "/home/ubuntu/data/candidates.bam"

println "comparing input bam: $params.input_bam with $params.truth_vcf"
println "outdir: $params.outdir"

ref = file(params.reference_fasta)

input_depth = 30 // could calculate this from inputs, or check that the inputs are consistent depth
bams = Channel.fromPath(params.input_bam)
depths = Channel.from(5,10,20,30)
bam_depths = depths.combine(bams).map{ depth, bam -> tuple(depth,"${bam.baseName}",bam) } 

process downsample{
    input:
    set depth, basename, bam from bam_depths

    output:
    set depth, basename, file("*x.bam") into bam_sniffles, bam_pbsv

    script:
    '''
    frac=$( bc -l <<< "scale=16; !{depth}/!{input_depth}" )

    if [[ frac < 0.9 ]]; then 
        sambamba view -t 8 -f bam -s $frac -o "!{basename}.${depth}x.bam" "!{bam}"
    else
        ln -s "!{bam} "!{basename}_!{depth}x.bam"
    fi
    '''
}

/*
 * sniffles variant caller
 */
process sniffles {
    publishDir "${params.outdir}/sniffles", mode:'copy'

    input:
    set depth, basename, file(bam) from bam_sniffles

    output:
    set depth, basename, "sniffles", file("*.vcf.gz") into sniffles_vcf_ch

    script:
    '''
    sniffles -s 3 --skip_parameter_estimation -m !{bam} -v "!{basename}.!{depth}x.Sniffles_s3_ignoreParam.vcf"
    gzip "*.vcf"
    '''
}
/*
 * pbsv variant caller
 */
process pbsv {
    publishDir "${params.outdir}/pbsv", mode:'copy'

    input:
    set depth, basename, file(bam) from  bam_pbsv

    output:
    set depth, basename, "pbsv", file("*.svsig.gz") into pbsv_vcf_ch

    script:
    '''
    pbsv discover !{bam} !{basename}.!{depth}x.pbsv.svsig.gz
    pbsv call !{ref} !{basename}.!{depth}x.pbsv.svsig.gz !{basename}.!{depth}x.pbsv.vcf
    '''
}

/*
 * truevari compares one vcf with another to generate true positive, false positive and true negative rates
 */
process truevari {
    publishDir "${params.outdir}/truevari", mode:'copy'

    input:
    set depth, basename, caller, file(vcf) from pbsv_vcf_ch.mix(sniffles_vcf_ch)

    output:
    file '*' into truevari_out_channel

    script:
    '''
    echo !{params.truth_vcf} 
    touch "!{basename}.!{depth}.!{caller}.truvari"
    '''
}
