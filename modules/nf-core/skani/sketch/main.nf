process SKANI_SKETCH {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/skani:0.2.2--ha6fb395_2':
        'biocontainers/skani:0.2.2--ha6fb395_2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}")                 , emit: sketch_dir
    tuple val(meta), path("${prefix}/${fasta}.sketch") , emit: sketch
    tuple val(meta), path("${prefix}/markers.bin")     , emit: markers
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    skani \\
        sketch \\
            ${fasta} \\
            -o ${prefix} \\
            -t ${task.cpus} \\
            ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: \$(skani --version 2>&1 | sed 's/^.*skani //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${fasta}.sketch
    touch ${prefix}/markers.bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: \$(skani --version 2>&1 | sed 's/^.*skani //; s/ .*\$//')
    END_VERSIONS
    """
}
