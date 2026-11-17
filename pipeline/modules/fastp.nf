// Adapter and quality trimming
process FASTP {
    tag "${meta.id}"
    label 'process_medium'
    publishDir path: { "${params.outdir}/${meta.id}/qc/fastp" }, mode: 'copy', pattern: '*.{json,html}'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_R1.trim.fastq.gz"), path("${meta.id}_R2.trim.fastq.gz"), emit: reads
    path "${meta.id}.fastp.json", emit: json
    path "${meta.id}.fastp.html", emit: html

    script:
    """
    fastp \\
        -i ${r1} -I ${r2} \\
        -o ${meta.id}_R1.trim.fastq.gz -O ${meta.id}_R2.trim.fastq.gz \\
        --detect_adapter_for_pe \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html
    """
}
