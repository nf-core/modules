process PROTEINORTHO {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/proteinortho:6.3.0--h70414c8_0':
        'quay.io/biocontainers/proteinortho:6.3.0--h70414c8_0' }"

    input:
    tuple val(meta), path(fasta_files, stageAs: "?/*")

    output:
    tuple val(meta), path("${prefix}.proteinortho.tsv")  , emit: orthologgroups
    tuple val(meta), path("${prefix}.proteinortho-graph"), emit: orthologgraph
    tuple val(meta), path("${prefix}.blast-graph")       , emit: blastgraph
    tuple val("${task.process}"), val('proteinortho'), eval("proteinortho --version 2>&1")                            , topic: versions, emit: versions_proteinortho
    tuple val("${task.process}"), val('diamond')     , eval("diamond version 2>/dev/null | sed '1!d;s/^.*version //'"), topic: versions, emit: versions_diamond
    tuple val("${task.process}"), val('blast')       , eval("blastp -version 2>/dev/null | sed '1!d;s/^.*: //;s/+//'"), topic: versions, emit: versions_blast

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_sorted = fasta_files.sort { a, b -> a.name <=> b.name }.join(' ')
    """
    proteinortho \\
        ${args} \\
        -cpus=${task.cpus} \\
        -project=${prefix} \\
        ${fasta_sorted}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.proteinortho.tsv
    touch ${prefix}.proteinortho-graph
    touch ${prefix}.blast-graph
    """
}
