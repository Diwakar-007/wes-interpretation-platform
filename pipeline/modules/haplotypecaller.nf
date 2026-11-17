// Call variants on one shard of the exome, in gVCF mode (allows joint calling of trios later)
process HAPLOTYPECALLER {
    tag "${meta.id}:${intervals.baseName}"

    input:
    tuple val(meta), path(bam), path(bai), path(intervals)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("${meta.id}.${intervals.baseName}.g.vcf.gz"), path("${meta.id}.${intervals.baseName}.g.vcf.gz.tbi"), emit: gvcf

    script:
    def mem = task.memory.toGiga() - 1
    """
    gatk --java-options "-Xmx${mem}g" HaplotypeCaller \\
        -R ${fasta} -I ${bam} \\
        -L ${intervals} \\
        -ERC GVCF \\
        --native-pair-hmm-threads ${task.cpus} \\
        -O ${meta.id}.${intervals.baseName}.g.vcf.gz
    """
}
