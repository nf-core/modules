process WISECONDORX_PREDICT {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::wisecondorx=1.2.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wisecondorx:1.2.5--pyh5e36f6f_0':
        'biocontainers/wisecondorx:1.2.5--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(npz)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(blacklist)

    output:
    tuple val(meta), path("*.bed"), emit: bed, optional:true
    tuple val(meta), path("*.png"), emit: plot, optional:true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--bed --plot'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed = blacklist ? "--blacklist ${blacklist}" : ""

    def VERSION = '1.2.5' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    WisecondorX predict \\
        ${npz} \\
        ${reference} \\
        ${prefix} \\
        ${bed} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisecondorx: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: '--bed --plot'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.5' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def bed = args.contains("--bed") ? "touch ${prefix}.bed" : ""
    def plot = args.contains("--plot") ? "touch ${prefix}.png" : ""

    """
    ${bed}
    ${plot}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisecondorx: ${VERSION}
    END_VERSIONS
    """
}
