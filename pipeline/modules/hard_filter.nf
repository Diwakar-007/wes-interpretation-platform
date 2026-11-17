// GATK-recommended hard filters (VQSR needs many samples, so single exomes use hard filters),
// then split multi-allelic sites and left-align indels so annotation is consistent.
// Failing variants are kept but tagged in the FILTER column.
process HARD_FILTER {
    tag "${meta.id}"
    label 'process_low'
    publishDir path: { "${params.outdir}/${meta.id}/vcf" }, mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), path("${meta.id}.filtered.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk SelectVariants -V ${vcf} -select-type SNP -O snps.vcf.gz
    gatk SelectVariants -V ${vcf} -select-type INDEL -select-type MIXED -O indels.vcf.gz

    gatk VariantFiltration -V snps.vcf.gz -O snps.filt.vcf.gz \\
        -filter "QD < 2.0"              --filter-name "QD2" \\
        -filter "QUAL < 30.0"           --filter-name "QUAL30" \\
        -filter "SOR > 3.0"             --filter-name "SOR3" \\
        -filter "FS > 60.0"             --filter-name "FS60" \\
        -filter "MQ < 40.0"             --filter-name "MQ40" \\
        -filter "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \\
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

    gatk VariantFiltration -V indels.vcf.gz -O indels.filt.vcf.gz \\
        -filter "QD < 2.0"               --filter-name "QD2" \\
        -filter "QUAL < 30.0"            --filter-name "QUAL30" \\
        -filter "FS > 200.0"             --filter-name "FS200" \\
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"

    gatk MergeVcfs -I snps.filt.vcf.gz -I indels.filt.vcf.gz -O merged.filt.vcf.gz

    bcftools norm -m -any -f ${fasta} -Oz -o ${meta.id}.filtered.vcf.gz merged.filt.vcf.gz
    bcftools index -t ${meta.id}.filtered.vcf.gz
    """
}
