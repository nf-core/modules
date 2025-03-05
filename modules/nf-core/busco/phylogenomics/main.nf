process BUSCO_PHYLOGENOMICS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco_phylogenomics:20240919--pyhdfd78af_0':
        'biocontainers/busco_phylogenomics:20240919--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(busco)

    output:
    tuple val(meta), path("${prefix}"), emit: phylogenomics
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '20240919' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    BUSCO_phylogenomics.py \\
        -i ${busco} \\
        -o ${prefix} \\
        -t $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco_phylogenomics: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '20240919' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mkdir ${prefix}
    touch ${prefix}/${prefix}.tree

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco_phylogenomics: $VERSION
    END_VERSIONS
    """
}
