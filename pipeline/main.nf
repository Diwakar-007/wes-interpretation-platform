#!/usr/bin/env nextflow
/*
 * Germline WES pipeline: FASTQ -> filtered, normalised VCF
 *   PREPARE_INTERVALS : pad capture targets, split into shards
 *   PREPROCESS        : FastQC, fastp, bwa mem, MarkDuplicates, BQSR
 *   CALL_VARIANTS     : HaplotypeCaller (gVCF, scattered) -> GenotypeGVCFs -> hard filters -> bcftools norm
 *   VEP_ANNOTATE      : Ensembl VEP + ClinVar + AlphaMissense + REVEL (only if --vep_cache is set)
 *   QC                : mosdepth coverage + MultiQC report
 */

include { PREPARE_INTERVALS } from './modules/prepare_intervals'
include { VEP_ANNOTATE      } from './modules/vep'
include { MOSDEPTH          } from './modules/mosdepth'
include { MULTIQC           } from './modules/multiqc'
include { PREPROCESS        } from './subworkflows/preprocess'
include { CALL_VARIANTS     } from './subworkflows/call_variants'

workflow {
    main:
    if (!params.input)       { error "Missing --input (samplesheet CSV)" }
    if (!params.fasta)       { error "Missing --fasta (reference genome)" }
    if (!params.known_sites) { error "Missing --known_sites (BQSR VCFs)" }
    if (!params.intervals)   { error "Missing --intervals (capture target BED)" }

    // ---------- inputs ----------
    ch_reads = channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            tuple([id: row.sample],
                  file(row.fastq_1, checkIfExists: true),
                  file(row.fastq_2, checkIfExists: true))
        }

    // ---------- reference files (value channels: reused by every task) ----------
    ch_fasta     = channel.value(file(params.fasta, checkIfExists: true))
    ch_fai       = channel.value(file("${params.fasta}.fai", checkIfExists: true))
    ch_dict      = channel.value(file(params.fasta.replaceFirst(/\.(fa|fasta)$/, '.dict'), checkIfExists: true))
    ch_bwa_index = channel.value(files("${params.fasta}.{amb,ann,bwt,pac,sa,alt}"))
    ch_bed       = channel.value(file(params.intervals, checkIfExists: true))

    def known_vcfs = params.known_sites.tokenize(',')
    ch_known     = channel.value(known_vcfs.collect { f -> file(f, checkIfExists: true) })
    ch_known_tbi = channel.value(known_vcfs.collect { f -> file("${f}.tbi", checkIfExists: true) })

    // ---------- workflow ----------
    PREPARE_INTERVALS(ch_bed, ch_fasta, ch_fai, ch_dict)

    PREPROCESS(ch_reads, ch_fasta, ch_fai, ch_dict, ch_bwa_index,
               ch_known, ch_known_tbi, PREPARE_INTERVALS.out.padded)

    CALL_VARIANTS(PREPROCESS.out.bam, ch_fasta, ch_fai, ch_dict,
                  PREPARE_INTERVALS.out.scattered)

    MOSDEPTH(PREPROCESS.out.bam, ch_bed)

    ch_qc = PREPROCESS.out.qc.mix(CALL_VARIANTS.out.qc, MOSDEPTH.out.qc)

    if (params.vep_cache) {
        VEP_ANNOTATE(
            CALL_VARIANTS.out.vcf, ch_fasta, ch_fai,
            channel.value(file(params.vep_cache, checkIfExists: true)),
            channel.value(file(params.vep_plugins, checkIfExists: true)),
            channel.value(tuple(file(params.clinvar, checkIfExists: true), file("${params.clinvar}.tbi", checkIfExists: true))),
            channel.value(tuple(file(params.alphamissense, checkIfExists: true), file("${params.alphamissense}.tbi", checkIfExists: true))),
            channel.value(tuple(file(params.revel, checkIfExists: true), file("${params.revel}.tbi", checkIfExists: true)))
        )
        ch_qc = ch_qc.mix(VEP_ANNOTATE.out.report)
    }

    MULTIQC(ch_qc.collect())
}
