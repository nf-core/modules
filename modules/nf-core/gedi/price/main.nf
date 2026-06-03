process GEDI_PRICE {
    tag "$meta.id"
    label 'process_medium'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cd/cd008e5721759d5909909254c77ec449778e0fc7c669b7c926b68f0c9059f510/data' :
        'community.wave.seqera.io/library/gedi_price:2392624d5f803049' }"

    input:
    tuple val(meta), path(bams, stageAs: 'bams/*'), path(bais, stageAs: 'bams/*')
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("${prefix}.orfs.tsv")                                                 , emit: orfs_tsv
    tuple val(meta), path("${prefix}.orfs.cit")                                                 , emit: orfs_cit, optional: true
    tuple val(meta), path("${prefix}.orfs.cit.metadata.json")                                   , emit: orfs_metadata, optional: true
    tuple val(meta), path("${prefix}.codons.cit")                                               , emit: codons_cit, optional: true
    tuple val(meta), path("${prefix}.model")                                                    , emit: model, optional: true
    tuple val(meta), path("${prefix}.signal.tsv")                                               , emit: signal, optional: true
    tuple val(meta), path("${prefix}.param")                                                    , emit: param, optional: true
    tuple val("${task.process}"), val('gedi'), eval("gedi -e Version 2>&1 | sed -n 's/.*Gedi version \\([^ ]*\\).*/\\1/p' | head -n 1"), topic: versions, emit: versions_gedi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def oml = "${index}/${meta2.id ?: 'reference'}.oml"
    """
    ls -1 bams/*.bam > price_input.bamlist
    bamlist2cit -n ${task.cpus} -p price_input.bamlist

    gedi -e Price \\
        -reads price_input.bamlist.cit \\
        -genomic ${oml} \\
        -prefix ${prefix} \\
        -nthreads ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.orfs.tsv
    touch ${prefix}.orfs.cit
    touch ${prefix}.orfs.cit.metadata.json
    touch ${prefix}.codons.cit
    touch ${prefix}.model
    touch ${prefix}.signal.tsv
    touch ${prefix}.param
    """
}
