process ANARCII {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/python_pip_anarcii:4e5c3ffabd22d3fc':
    'community.wave.seqera.io/library/python_pip_anarcii:4e5c3ffabd22d3fc' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.csv"), emit: anarcii
    tuple val("${task.process}"), val('anarcii'), eval('anarcii --version | sed "s/^anarcii //"'), topic: versions, emit: versions_anarcii

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    anarcii \\
        $args\\
        -o ${prefix}.csv \\
        ${fasta}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    touch ${prefix}.csv
    """
}
