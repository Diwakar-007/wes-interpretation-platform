// Coverage over capture targets: mean depth and % of bases at >=10x/20x/30x/50x
process MOSDEPTH {
    tag "${meta.id}"
    label 'process_low'
    publishDir path: { "${params.outdir}/${meta.id}/qc/mosdepth" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path bed

    output:
    path "${meta.id}.*", emit: qc

    script:
    """
    mosdepth --threads ${task.cpus} --by ${bed} --thresholds 10,20,30,50 \\
        --no-per-base --fast-mode ${meta.id} ${bam}
    """
}
