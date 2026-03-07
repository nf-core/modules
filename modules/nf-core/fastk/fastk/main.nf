process FASTK_FASTK {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/02/02c05b2ec421debc83883ef9a211291e3220546c12f7c54cb78e66209cb2797d/data' :
        'community.wave.seqera.io/library/fastk:1.2--4bc70c6cd0d420bd' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.hist")                      , emit: hist
    tuple val(meta), path("*.ktab*", hidden: true)       , emit: ktab, optional: true
    tuple val(meta), path("*.{prof,pidx}*", hidden: true), emit: prof, optional: true
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('fastk'), val('1.2'), emit: versions_fastk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    """
    FastK \\
        $args \\
        -T$task.cpus \\
        -M${task.memory.toGiga()} \\
        -N${prefix} \\
        $reads

    find . -name '*.ktab*' -exec chmod a+r {} \\;
    """

    stub:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def touch_ktab = args.contains('-t') ? "touch ${prefix}.ktab .${prefix}.ktab.1" : ''
    def touch_prof = args.contains('-p') ? "touch ${prefix}.prof .${prefix}.pidx.1" : ''
    """
    touch ${prefix}.hist
    $touch_ktab
    $touch_prof

    echo \\
    "FastK \\
        $args \\
        -T$task.cpus \\
        -M${task.memory.toGiga()} \\
        -N${prefix}_fk \\
        $reads"
    """
}
