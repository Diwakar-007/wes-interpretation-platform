// Turn the gVCF (all sites) into a VCF of variant sites with genotypes
process GENOTYPE_GVCFS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(gvcf), path(tbi)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("${meta.id}.raw.vcf.gz"), path("${meta.id}.raw.vcf.gz.tbi"), emit: vcf

    script:
    def mem = task.memory.toGiga() - 1
    """
    gatk --java-options "-Xmx${mem}g" GenotypeGVCFs \\
        -R ${fasta} -V ${gvcf} \\
        -O ${meta.id}.raw.vcf.gz
    """
}
