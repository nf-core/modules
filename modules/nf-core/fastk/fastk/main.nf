process FASTK_FASTK {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9f/9f0bee9bfacd05665a9b1a11dd087dbf1be41ac3e640931c38c914a2390642cf/data' :
        'community.wave.seqera.io/library/fastk_merquryfk_r-cowplot_r-ggplot2_r-viridis:f9994edc2270683c' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.hist")                      , emit: hist
    tuple val(meta), path("*.ktab*", hidden: true)       , emit: ktab, optional: true
    tuple val(meta), path("*.{prof,pidx}*", hidden: true), emit: prof, optional: true
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FASTK_VERSION = '1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    FastK \\
        $args \\
        -T$task.cpus \\
        -M${task.memory.toGiga()} \\
        -N${prefix}_fk \\
        $reads

    find . -name '*.ktab*' \\
        | xargs chmod a+r

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FASTK_VERSION = '1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def touch_ktab = args.contains('-t') ? "touch ${prefix}_fk.ktab .${prefix}_fk.ktab.1" : ''
    def touch_prof = args.contains('-p') ? "touch ${prefix}_fk.prof .${prefix}_fk.pidx.1" : ''
    """
    touch ${prefix}_fk.hist
    $touch_ktab
    $touch_prof

    echo \\
    "FastK \\
        $args \\
        -T$task.cpus \\
        -M${task.memory.toGiga()} \\
        -N${prefix}_fk \\
        $reads"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
    END_VERSIONS
    """
}
