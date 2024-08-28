process CNVPYTOR_IMPORTREADDEPTH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvpytor:1.3.1--pyhdfd78af_1':
        'biocontainers/cnvpytor:1.3.1--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(input_file), path(index)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.pytor")	, emit: pytor
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "-T ${fasta}" : ''
    """
    cnvpytor \\
        -root ${prefix}.pytor \\
        -rd $input_file \\
        $args \\
        $reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(cnvpytor --version | sed -n 's/.*CNVpytor \\(.*\\)/\\1/p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pytor

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(cnvpytor --version | sed -n 's/.*CNVpytor \\(.*\\)/\\1/p')
    END_VERSIONS
    """
}
