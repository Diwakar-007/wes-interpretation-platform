// Variant summary statistics (counts, Ts/Tv ratio, indel sizes) for MultiQC
process BCFTOOLS_STATS {
    tag "${meta.id}"
    label 'process_low'
    publishDir path: { "${params.outdir}/${meta.id}/qc/bcftools" }, mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path fasta

    output:
    path "${meta.id}.bcftools_stats.txt", emit: stats

    script:
    """
    bcftools stats -F ${fasta} -f PASS ${vcf} > ${meta.id}.bcftools_stats.txt
    """
}
