process FOLDCOMP_DECOMPRESS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/foldcomp:1.0.0--h7f5d12c_0':
        'biocontainers/foldcomp:1.0.0--h7f5d12c_0' }"

    input:
    tuple val(meta), path(fcz)

    output:
    tuple val(meta), path("{*pdb,*.cif}"), emit: pdb
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    foldcomp \\
        $args \\
        decompress \\
        -t ${task.cpus} \\
        ${fcz}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldcomp: \$(foldcomp --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pdb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldcomp: \$(foldcomp --version)
    END_VERSIONS
    """
}
