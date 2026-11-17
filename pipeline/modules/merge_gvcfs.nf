// Gather the per-shard gVCFs of one sample into a single gVCF
process MERGE_GVCFS {
    tag "${meta.id}"
    label 'process_low'
    publishDir path: { "${params.outdir}/${meta.id}/vcf" }, mode: 'copy'

    input:
    tuple val(meta), path(gvcfs), path(tbis)
    path dict

    output:
    tuple val(meta), path("${meta.id}.g.vcf.gz"), path("${meta.id}.g.vcf.gz.tbi"), emit: gvcf

    script:
    def inputs = [gvcfs].flatten().collect { g -> "-I ${g}" }.join(' ')
    """
    gatk --java-options "-Xmx1g" MergeVcfs \\
        ${inputs} \\
        -D ${dict} \\
        -O ${meta.id}.g.vcf.gz
    """
}
