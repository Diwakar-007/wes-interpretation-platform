// Pad capture targets once, then split them into shards for parallel calling.
// Padding BEFORE splitting avoids overlapping shards (which would duplicate variants).
process PREPARE_INTERVALS {
    label 'process_low'

    input:
    path bed
    path fasta
    path fai
    path dict

    output:
    path "padded.interval_list",   emit: padded
    path "scatter/*.interval_list", emit: scattered

    script:
    """
    gatk --java-options "-Xmx1g" PreprocessIntervals \\
        -R ${fasta} -L ${bed} \\
        --padding ${params.interval_padding} \\
        --bin-length 0 \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -O padded.interval_list

    gatk --java-options "-Xmx1g" SplitIntervals \\
        -R ${fasta} -L padded.interval_list \\
        --scatter-count ${params.scatter_count} \\
        --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \\
        -O scatter
    """
}
