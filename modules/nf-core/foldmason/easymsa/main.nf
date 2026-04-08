process FOLDMASON_EASYMSA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/29/299b6d1a8a81632118e9c0cd110f9d89aae58b5a0132543a7fedbecb656da570/data':
            'community.wave.seqera.io/library/foldmason_pigz:289a5a7a978299df' }"


    input:
    tuple val(meta) , path(pdbs)
    tuple val(meta2), path(tree)
    val(compress)

    output:
    tuple val(meta), path("${prefix}_3di.fa${compress ? '.gz' : ''}"), emit: msa_3di
    tuple val(meta), path("${prefix}_aa.fa${compress ? '.gz' : ''}") , emit: msa_aa
    tuple val("${task.process}"), val('foldmason'), eval('foldmason version'), emit: versions_foldmason, topic: versions
    tuple val("${task.process}"), val('pigz'), eval('pigz --version 2>&1 | sed "s/^.*pigz[[:space:]]*//"'), emit: versions_pigz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def options_tree = tree ? "--guide-tree $tree" : ""
    """
    foldmason easy-msa \\
        ${pdbs} \\
        ${prefix} \\
        tmp \\
        ${options_tree} \\
        $args \\
        --threads $task.cpus

    if ${compress}; then
        pigz -p ${task.cpus} ${prefix}_3di.fa
        pigz -p ${task.cpus} ${prefix}_aa.fa
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ""  ${compress ? '| gzip' : ''} > ${prefix}_3di.fa${compress ? '.gz' : ''}
    echo ""  ${compress ? '| gzip' : ''} > ${prefix}_aa.fa${compress ? '.gz' : ''}
    """
}
