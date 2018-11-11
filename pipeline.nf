/*
 * pipeline input parameters
 */
params.input_bams = "/home/ubuntu/data/HG002.10kb.Q20.GRCh38.pbmm2.bam"
params.outdir = "/home/ubuntu/pipeline-output"
params.reference_fasta_path = "/home/ubuntu/data/reference/hs37d5.fa"
params.bed = "/home/ubuntu/data/candidates.bam"

println "comparing input bams: $params.input_bams with $params.truth_vcf"
println "outdir: $params.outdir"

ref = file(params.reference_fasta_path)

input_depth = 30 // could calculate this from inputs, or check that the inputs are consistent depth
bams = Channel.fromPath(params.input_bams)
depths = Channel.from(5,10,20,30)
bam_depths = depths.combine(bams).map{ depth, bam -> tuple(depth,"${bam.baseName}",file(bam)) } 

process downsample{
    input:
    set depth, basename, bam from bam_depths

    output:
    set depth, basename, file("*x.bam") into bam_sniffles, bam_pbsv

    shell:
    '''
    frac=$( bc -l <<< "scale=16; !{depth}/!{input_depth}" )

    if [[ $frac < 0.9 ]]; then 
        sambamba view -t 8 -f bam -s $frac -o "!{basename}.!{depth}x.bam" "!{bam}"
    else
        ln -s "!{bam}" "!{basename}_!{depth}x.bam"
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
    set depth, basename, val("sniffles"), file("*.vcf.gz"),file("*.vcf.gz.tbi") into sniffles_vcf_ch

    shell:
    '''
    sniffles -s 3 --skip_parameter_estimation -m !{bam} -v "!{basename}.!{depth}x.Sniffles_s3_ignoreParam.vcf"
    bgzip *.vcf
    tabix -p vcf *.vcf.gz
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
    set depth, basename, val("pbsv"), file("*.pbsv.vcf.gz"), file("*.pbsv.vcf.gz.tbi") into pbsv_vcf_ch

    shell:
    '''
    pbsv discover !{bam} !{basename}.!{depth}x.pbsv.svsig.gz
    pbsv call !{ref} !{basename}.!{depth}x.pbsv.svsig.gz !{basename}.!{depth}x.pbsv.vcf
    bgzip *.vcf
    tabix -p vcf *.vcf.gz
    '''
}

/*
 * truevari compares one vcf with another to generate true positive, false positive and true negative rates
 */
process truvari {
     publishDir "${params.outdir}/truvari", mode:'copy'

     input:
     set depth, basename, caller, file(vcf), file(tbi) from pbsv_vcf_ch.mix(sniffles_vcf_ch)

     output:
     file 'out/*giab_report.txt' into truvari_giab_ch
     file 'out/*tp-base.vcf' into truvari_tp_base_ch

     shell:
     '''
     truvari --passonly -b !{params.truth_vcf} -c !{vcf} -o out -f !{ref} --giabreport
     mv out/giab_report.txt out/!{basename}.!{depth}.!{caller}.giab_report.txt
     mv out/tp-base.vcf out/!{basename}.!{depth}.!{caller}.tp-base.vcf
     '''
}

/*
 * coalesces truvari true positive vcfs to form a big comparison set for visualization
 */
process combine_truvari_tp {
    publishDir "${params.outdir}/truvari", mode: 'copy'
    input:
    file ('*') from truvari_tp_base_ch.flatten().toList()

    output:
    file "upset_input.txt" into upset_input_ch

    shell:
    '''
    cat * > upset_input.txt
    '''
}

/*
* coalesces truvari giab files into a single report
*/
process combine_truvari_giab {
    publishDir "${params.outdir}/truvari", mode: 'copy'
    input:
    file ('*') from truvari_giab_ch.flatten().toList()

    output:
    file "all_truvari_comparisons.tsv" into combined_giab_ch

    shell:
    '''
    cat * > all_truvari_comparisons.tsv
    '''
}
