// Annotate variants with Ensembl VEP (offline cache) + ClinVar + AlphaMissense + REVEL.
// All transcripts are kept; one consequence per allele+gene is flagged PICK=1, preferring
// the MANE Select transcript (the transcript clinical reports use).
process VEP_ANNOTATE {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/${meta.id}/vcf" }, mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path fasta
    path fai
    path vep_cache
    path vep_plugins
    tuple path(clinvar), path(clinvar_tbi)
    tuple path(alphamissense), path(alphamissense_tbi)
    tuple path(revel), path(revel_tbi)

    output:
    tuple val(meta), path("${meta.id}.vep.vcf.gz"), path("${meta.id}.vep.vcf.gz.tbi"), emit: vcf
    path "${meta.id}.vep_summary.html", emit: report

    script:
    """
    vep \\
        --input_file ${vcf} --format vcf \\
        --output_file ${meta.id}.vep.vcf.gz --vcf --compress_output bgzip \\
        --stats_file ${meta.id}.vep_summary.html \\
        --offline --cache --dir_cache ${vep_cache} --cache_version ${params.vep_cache_version} \\
        --assembly GRCh38 --fasta ${fasta} \\
        --fork ${task.cpus} \\
        --hgvs --hgvsg --symbol --biotype --canonical --mane --numbers --protein --ccds \\
        --variant_class --domains \\
        --check_existing --af_gnomade --af_gnomadg --max_af \\
        --flag_pick_allele_gene \\
        --pick_order mane_select,canonical,appris,tsl,biotype,ccds,rank,length \\
        --custom file=${clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNSIGCONF \\
        --dir_plugins ${vep_plugins} \\
        --plugin AlphaMissense,file=${alphamissense} \\
        --plugin REVEL,file=${revel}

    tabix -p vcf ${meta.id}.vep.vcf.gz
    """
}
