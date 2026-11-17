// Align to the reference and coordinate-sort.
// -K fixes the batch size so results are identical regardless of thread count.
// -Y uses soft-clipping for supplementary alignments (GATK best practice).
process BWA_MEM {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(r1), path(r2)
    path fasta
    path fai
    path bwa_index

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), emit: bam

    script:
    def read_group = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}"
    """
    bwa mem -K 100000000 -Y -t ${task.cpus} -R '${read_group}' ${fasta} ${r1} ${r2} \\
        | samtools sort -@ 2 -m 768M -o ${meta.id}.sorted.bam -
    """
}
