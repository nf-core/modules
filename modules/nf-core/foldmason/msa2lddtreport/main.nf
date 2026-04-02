process FOLDMASON_MSA2LDDTREPORT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dc0b27a73277c851f0a1df09fc5fa4f0ff80832e6779f143420a7e2321ac62b4/data' :
            'community.wave.seqera.io/library/foldmason:4.dd3c235--375019dcc3c5fe6f' }"
    input:
    tuple val(meta) , path(msa)
    tuple val(meta2), path(db)
    tuple val(meta3), path(pdbs)
    tuple val(meta4), path(tree)

    output:
    tuple val(meta), path("${prefix}.html"), emit: html
    tuple val("${task.process}"), val('foldmason'), eval('foldmason version'), emit: versions_foldmason, topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def options_tree = tree ? "--guide-tree $tree" : ""
    """
    foldmason msa2lddtreport \\
        ${meta.id} \\
        ${msa} \\
        ${prefix}.html \\
        $args \\
        ${options_tree} \\
        --threads $task.cpus
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html
    """
}
