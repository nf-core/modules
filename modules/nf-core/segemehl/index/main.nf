process SEGEMEHL_INDEX {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/segemehl:0.3.4--hc2ea5fd_5':
        'quay.io/biocontainers/segemehl:0.3.4--hc2ea5fd_5' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.idx"), emit: index
    tuple val("${task.process}"), val('segemehl'), eval("segemehl.x 2>&1 | sed '/^  [0-9]\\.[0-9\\.]*/!d;s/^  //;s/ .*//'"), topic: versions, emit: versions_segemehl

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    segemehl.x \\
        -t ${task.cpus} \\
        -d ${fasta} \\
        -x ${prefix}.idx \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    touch ${prefix}.idx
    """
}
