process FGBIO_COLLECTDUPLEXSEQMETRICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d6/d65e7c77d89e7e443384d17a94ffe31fe988b05bc7d695f2a75beaf502721925/data':
        'community.wave.seqera.io/library/fgbio_r-ggplot2:cf2b9a5308d77b67' }"

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
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals = interval_list ? "--intervals ${bed}" : ""
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
    END_VERSIONS
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
    END_VERSIONS
    """
}
