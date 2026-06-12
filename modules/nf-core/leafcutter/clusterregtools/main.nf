process LEAFCUTTER_CLUSTERREGTOOLS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/leafcutter:2.0.3-pyhd8ed1ab_0':
        'quay.io/biocontainers/leafcutter:2.0.3-pyhd8ed1ab_0' }"

    input:
    tuple val(meta), path(juncfiles)

    output:
    tuple val(meta), path("*_pooled")                      , emit: pooled
    tuple val(meta), path("*_refined")                     , emit: refined
    tuple val(meta), path("*_sortedlibs")                  , emit: sortedlibs
    tuple val(meta), path("*_perind*.counts.gz")           , emit: counts
    tuple val(meta), path("*_perind_numers*.counts.gz")    , emit: numers
    tuple val("${task.process}"), val('leafcutter'), eval("leafcutter-cluster --version 2>&1 | grep -oP 'leafcutter.*' || echo 'leafcutter 2.0.2'"), topic: versions, emit: versions_leafcutter

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    leafcutter-cluster \\
        -j $juncfiles \\
        -r ./ \\
        -o $prefix \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def const_suffix = task.ext.args?.contains("--includeconst") || task.ext.args?.contains("-C") ? ".const" : ""
    """
    touch ${prefix}_pooled
    touch ${prefix}_refined
    touch ${prefix}_sortedlibs
    echo "" | gzip > ${prefix}_perind${const_suffix}.counts.gz
    echo "" | gzip > ${prefix}_perind_numers${const_suffix}.counts.gz
    """
}
