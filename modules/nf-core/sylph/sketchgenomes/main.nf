process SYLPH_SKETCHGENOMES {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sylph:0.9.0--ha6fb395_0'
        : 'quay.io/biocontainers/sylph:0.9.0--ha6fb395_0'}"

    input:
    tuple val(meta), path(fasta, stageAs: 'genomes/')

    output:
    tuple val(meta), path('*.syldb'), emit: syldb
    tuple val("${task.process}"), val('sylph'), eval('sylph -V | sed "s/sylph //g"'), topic: versions, emit: versions_sylph

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    find -L genomes/ -type f > genomes.txt

    sylph sketch \\
        -t ${task.cpus} \\
        ${args} \\
        --gl genomes.txt \\
        -o ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    echo "${args}"
    touch ${prefix}.syldb
    """
}
