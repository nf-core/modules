process EMBOSS_REVSEQ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/emboss:6.6.0--h86d058a_5':
        'quay.io/biocontainers/emboss:6.6.0--h86d058a_5' }"

    input:
    tuple val(meta), path(sequences)

    output:
    tuple val(meta), path("*.${sequences.name - ~/.*\./}"), emit: revseq
    tuple val("${task.process}"), val("emboss"), eval("revseq -version 2>&1 | sed 's/EMBOSS://'"), topic: versions, emit: versions_emboss

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = sequences.name - ~/.*\./
    def outfile = "${prefix}.rev.${suffix}"
    """
    revseq \\
        $args \\
        $sequences \\
        $outfile
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = sequences.name - ~/.*\./
    def outfile = "${prefix}.rev.${suffix}"
    """
    touch ${outfile}
    """
}
