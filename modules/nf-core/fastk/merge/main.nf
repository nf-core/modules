process FASTK_MERGE {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9f/9f0bee9bfacd05665a9b1a11dd087dbf1be41ac3e640931c38c914a2390642cf/data' :
        'community.wave.seqera.io/library/fastk_merquryfk_r-cowplot_r-ggplot2_r-viridis:f9994edc2270683c' }"

    input:
    tuple val(meta), path(hist), path(ktab), path(prof)

    output:
    tuple val(meta), path("*.hist")                      , emit: hist
    tuple val(meta), path("*.ktab*", hidden: true)       , emit: ktab, optional: true
    tuple val(meta), path("*.{prof,pidx}*", hidden: true), emit: prof, optional: true
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('fastk'), eval('echo 1.1'), emit: versions_fastk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    """
    Fastmerge \\
        $args \\
        -T$task.cpus \\
        ${prefix} \\
        $hist
    """

    stub:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def touch_hist    = args.contains('-h') ? "touch ${prefix}.hist"                   : ''
    def touch_ktab    = args.contains('-t') ? "touch ${prefix}.ktab .${prefix}.ktab.1" : ''
    """
    ${touch_hist}
    ${touch_ktab}
    """
}
