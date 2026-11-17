// BAM -> scattered gVCFs -> merged gVCF -> genotyped VCF -> hard-filtered, normalised VCF
include { HAPLOTYPECALLER } from '../modules/haplotypecaller'
include { MERGE_GVCFS     } from '../modules/merge_gvcfs'
include { GENOTYPE_GVCFS  } from '../modules/genotype_gvcfs'
include { HARD_FILTER     } from '../modules/hard_filter'
include { BCFTOOLS_STATS  } from '../modules/bcftools_stats'

workflow CALL_VARIANTS {
    take:
    bam          // tuple(meta, bam, bai)
    fasta
    fai
    dict
    scattered    // list of interval_list shards

    main:
    // one HaplotypeCaller task per (sample, shard)
    ch_hc_input = bam.combine(scattered.flatten())
    HAPLOTYPECALLER(ch_hc_input, fasta, fai, dict)

    // gather the shards of each sample back together
    MERGE_GVCFS(HAPLOTYPECALLER.out.gvcf.groupTuple(), dict)
    GENOTYPE_GVCFS(MERGE_GVCFS.out.gvcf, fasta, fai, dict)
    HARD_FILTER(GENOTYPE_GVCFS.out.vcf, fasta, fai, dict)
    BCFTOOLS_STATS(HARD_FILTER.out.vcf, fasta)

    emit:
    vcf = HARD_FILTER.out.vcf    // tuple(meta, vcf.gz, tbi)
    qc  = BCFTOOLS_STATS.out.stats
}
