// Base Quality Score Recalibration: model systematic sequencer errors using known
// variant sites (which are masked), then rewrite base qualities.
process BQSR {
    tag "${meta.id}"
    label 'process_medium'
    publishDir path: { "${params.outdir}/${meta.id}/bam" }, mode: 'copy', pattern: '*.{bam,bai}'
    publishDir path: { "${params.outdir}/${meta.id}/qc/bqsr" }, mode: 'copy', pattern: '*.table'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path dict
    path known_sites
    path known_sites_tbi
    path intervals

    output:
    tuple val(meta), path("${meta.id}.recal.bam"), path("${meta.id}.recal.bai"), emit: bam
    path "${meta.id}.recal.table", emit: table

    script:
    def mem = task.memory.toGiga() - 1
    def known = [known_sites].flatten().collect { vcf -> "--known-sites ${vcf}" }.join(' ')
    """
    gatk --java-options "-Xmx${mem}g" BaseRecalibrator \\
        -R ${fasta} -I ${bam} ${known} \\
        -L ${intervals} \\
        -O ${meta.id}.recal.table

    gatk --java-options "-Xmx${mem}g" ApplyBQSR \\
        -R ${fasta} -I ${bam} \\
        --bqsr-recal-file ${meta.id}.recal.table \\
        -O ${meta.id}.recal.bam
    """
}
