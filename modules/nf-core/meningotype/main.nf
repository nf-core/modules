process MENINGOTYPE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meningotype:0.8.5--pyhdfd78af_0' :
        'quay.io/biocontainers/meningotype:0.8.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val("meningotype"), eval("meningotype --version 2>&1 | sed 's/^.*meningotype v//'"), emit: versions_meningotype, topic: versions
    tuple val("${task.process}"), val("biopython"), eval("pip show biopython | sed -n 's/Version: //p'"), emit: versions_biopython, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    meningotype \\
        ${args} \\
        ${fasta} \\
        > ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
