process PASTY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pasty:1.0.0--hdfd78af_0':
        'quay.io/biocontainers/pasty:1.0.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.tsv")        , emit: tsv
    tuple val(meta), path("${prefix}.blastn.tsv") , emit: blast
    tuple val(meta), path("${prefix}.details.tsv"), emit: details
    tuple val("${task.process}"), val('pasty'), eval("pasty --version 2>&1 | sed 's/^.*pasty, version //;'"), topic: versions, emit: versions_pasty

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    pasty \\
        ${args} \\
        --prefix ${prefix} \\
        --assembly ${fasta}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.blastn.tsv
    touch ${prefix}.details.tsv
    """
}
