process FOLDMASON_EASYMSA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/06/067e6389ab95497b753ba1deabaa6acbce25b99c8cfcf39c06d5c1af42fd7751/data':
            'community.wave.seqera.io/library/foldmason_pigz:88809eb5649534b0' }"


    input:
    tuple val(meta) , path(pdbs)
    tuple val(meta2), path(tree)
    val(compress)

    output:
    tuple val(meta), path("${prefix}_3di.fa${compress ? '.gz' : ''}"), emit: msa_3di
    tuple val(meta), path("${prefix}_aa.fa${compress ? '.gz' : ''}") , emit: msa_aa
    path "versions.yml"                                              , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldmason: \$(foldmason | grep "foldmason Version:" | cut -d":" -f 2 | awk '{\$1=\$1;print}')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ""  ${compress ? '| gzip' : ''} > ${prefix}_3di.fa${compress ? '.gz' : ''}
    echo ""  ${compress ? '| gzip' : ''} > ${prefix}_aa.fa${compress ? '.gz' : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldmason: \$(foldmason | grep "foldmason Version:" | cut -d":" -f 2 | awk '{\$1=\$1;print}')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
}
