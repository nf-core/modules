process FGBIO_COLLECTDUPLEXSEQMETRICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe9479adc5e6e0a1c125d346fdfa0dd313834249e9c55c40e8d44ec3a48c6559/data' :
        'community.wave.seqera.io/library/fgbio:3.1.1--6c9a88faf1d62b6c' }"
    input:
    tuple val(meta), path(grouped_bam)
    path interval_list

    output:
    tuple val(meta), path("**.family_sizes.txt")        , emit: family_sizes
    tuple val(meta), path("**.duplex_family_sizes.txt") , emit: duplex_family_sizes
    tuple val(meta), path("**.duplex_yield_metrics.txt"), emit: duplex_yield_metrics
    tuple val(meta), path("**.umi_counts.txt")          , emit: umi_counts
    tuple val(meta), path("**.duplex_qc.pdf")           , emit: duplex_qc
    tuple val(meta), path("**.duplex_umi_counts.txt")   , emit: duplex_umi_counts, optional: true
    tuple val("${task.process}"), val('fgbio'), eval('fgbio --version 2>&1 | tr -d "[:cntrl:]" | sed -e "s/^.*Version: //;s/\\[.*$//"'), topic: versions, emit: versions_fgbio
    tuple val("${task.process}"), val('ggplot2'), eval('Rscript -e "cat(as.character(packageVersion(\'ggplot2\')))"'), topic: versions, emit: versions_ggplot2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals = interval_list ? "--intervals ${interval_list}" : ""
    def mem_gb = 8

    if (!task.memory) {
        log.info '[fgbio CollectDuplexSeqMetrics] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else if (mem_gb > task.memory.giga) {
        if (task.memory.giga < 2) {
            mem_gb = 1
        } else {
            mem_gb = task.memory.giga - 1
        }
    }

    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --async-io=true \\
        --compression=1 \\
        CollectDuplexSeqMetrics \\
        --input $grouped_bam \\
        --output ${prefix} \\
        $intervals \\
        $args
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def touch_duplex_umi = args.contains("--duplex-umi-counts") || args.contains("-u") ? "touch ${prefix}.duplex_umi_counts.txt" : ""

    """
    touch ${prefix}.family_sizes.txt
    touch ${prefix}.duplex_family_sizes.txt
    touch ${prefix}.duplex_yield_metrics.txt
    touch ${prefix}.umi_counts.txt
    touch ${prefix}.duplex_qc.pdf
    $touch_duplex_umi
    """
}
