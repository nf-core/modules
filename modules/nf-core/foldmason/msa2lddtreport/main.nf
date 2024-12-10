process FOLDMASON_MSA2LDDTREPORT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a8/a88d162c3f39a1518d48c3faec235e6fcde750586da868b62fc5f0a08a89aa9d/data' :
            'community.wave.seqera.io/library/foldmason:2.7bd21ed--e7f739473ad6578d' }"
    input:
    tuple val(meta)  , path(msa)
    tuple val(meta2) , path(db)
    tuple val(meta3) , path(pdbs)
    tuple val(meta4) , path(tree)

    output:
    tuple val(meta), path("${prefix}.html"), emit: html
    path "versions.yml"                    , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldmason: \$(foldmason | grep "foldmason Version:" | cut -d":" -f 2 | awk '{\$1=\$1;print}')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldmason: \$(foldmason | grep "foldmason Version:" | cut -d":" -f 2 | awk '{\$1=\$1;print}')
    END_VERSIONS
    """
}
