// Raw read quality report
process FASTQC {
    tag "${meta.id}"
    label 'process_low'
    publishDir path: { "${params.outdir}/${meta.id}/qc/fastqc" }, mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    path "*_fastqc.zip",  emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${r1} ${r2}
    """
}
