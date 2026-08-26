process SEQUENZAUTILS_GCWIGGLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py38h6ed170a_2' :
        'quay.io/biocontainers/sequenza-utils:3.0.0--py38h6ed170a_2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.wig.gz"), emit: wig
    tuple val("${task.process}"), val("sequenzautils"), eval("sequenza-utils 2>&1 | sed '/version/!d;s/.*version //;s/ .*//'"), topic: versions, emit: versions_sequenzautils

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sequenza-utils \\
        gc_wiggle \\
        ${args} \\
        --fasta ${fasta} \\
        -o ${prefix}.wig.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.wig.gz
    """
}
