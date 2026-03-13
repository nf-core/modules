process MERQURYFK_PLOIDYPLOT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7d113fd30596756ca41c1651aecff1dfde15d2f08dc34996e5b89c7eabd50ab1/data' :
        'community.wave.seqera.io/library/fastk_merquryfk_r-argparse_r-cowplot_pruned:985041255d013130' }"

    input:
    tuple val(meta), path(fastk_hist), path(fastk_ktab)

    output:
    tuple val(meta), path("*.{png,pdf}"), emit: images
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('merquryfk'), eval('echo 1.1.1'), emit: versions_merquryfk, topic: versions
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('fastk'), eval('echo 1.1'), emit: versions_fastk, topic: versions
    tuple val("${task.process}"), val('R'), eval('R --version | sed "1!d; s/.*version //; s/ .*//"'), emit: versions_r, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    """
    PloidyPlot \\
        $args \\
        -T$task.cpus \\
        ${fastk_ktab.find { path -> path.toString().endsWith(".ktab") }} \\
        ${prefix}
    """

    stub:
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    def outfmt          = args.contains('-pdf') ? "pdf" : "png"
    """
    touch ${prefix}.fi.${outfmt}
    touch ${prefix}.ln.${outfmt}
    touch ${prefix}.st.${outfmt}
    """
}
