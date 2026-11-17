// FASTQ -> QC -> trimmed reads -> aligned, deduplicated, recalibrated BAM
include { FASTQC         } from '../modules/fastqc'
include { FASTP          } from '../modules/fastp'
include { BWA_MEM        } from '../modules/bwa_mem'
include { MARKDUPLICATES } from '../modules/markduplicates'
include { BQSR           } from '../modules/bqsr'

workflow PREPROCESS {
    take:
    reads            // tuple(meta, r1, r2)
    fasta
    fai
    dict
    bwa_index
    known_sites
    known_sites_tbi
    intervals        // padded interval_list

    main:
    FASTQC(reads)
    FASTP(reads)
    BWA_MEM(FASTP.out.reads, fasta, fai, bwa_index)
    MARKDUPLICATES(BWA_MEM.out.bam)
    BQSR(MARKDUPLICATES.out.bam, fasta, fai, dict, known_sites, known_sites_tbi, intervals)

    qc = FASTQC.out.zip.mix(FASTP.out.json, MARKDUPLICATES.out.metrics, BQSR.out.table)

    emit:
    bam = BQSR.out.bam    // tuple(meta, bam, bai)
    qc  = qc
}
