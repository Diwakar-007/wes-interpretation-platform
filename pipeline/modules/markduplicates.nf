// Flag PCR/optical duplicates so the variant caller ignores them
process MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_medium'
    publishDir path: { "${params.outdir}/${meta.id}/qc/markduplicates" }, mode: 'copy', pattern: '*.txt'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.md.bam"), path("${meta.id}.md.bai"), emit: bam
    path "${meta.id}.md_metrics.txt", emit: metrics

    script:
    def mem = task.memory.toGiga() - 1
    """
    gatk --java-options "-Xmx${mem}g" MarkDuplicates \\
        -I ${bam} \\
        -O ${meta.id}.md.bam \\
        -M ${meta.id}.md_metrics.txt \\
        --CREATE_INDEX true \\
        --TMP_DIR .
    """
}
