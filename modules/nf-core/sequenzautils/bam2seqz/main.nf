process SEQUENZAUTILS_BAM2SEQZ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py38h6ed170a_2' :
        'quay.io/biocontainers/sequenza-utils:3.0.0--py38h6ed170a_2' }"

    input:
    tuple val(meta), path(normalbam), path(tumourbam)
    path fasta
    path wigfile

    output:
    tuple val(meta), path("*.gz"), emit: seqz
    tuple val("${task.process}"), val("sequenzautils"), eval("sequenza-utils 2>&1 | sed '/version/!d;s/.*version //;s/ .*//'"), topic: versions, emit: versions_sequenzautils

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sequenza-utils \\
        bam2seqz \\
        ${args} \\
        -n ${normalbam} \\
        -t ${tumourbam} \\
        --fasta ${fasta} \\
        -gc ${wigfile} \\
        -o ${prefix}.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.gz
    """
}
