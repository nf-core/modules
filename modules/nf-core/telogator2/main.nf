process TELOGATOR2 {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ab/ab4d9d463b2866006f8cbca9fbe6f978b1803e41f2a97d9f4d3c14ff6d97822f/data'
        : 'community.wave.seqera.io/library/telogator2:2.2.3--01b2748e09721f3b' }"

    input:
    tuple val(meta), path(reads), path(reads_index)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("${prefix}/tlens_by_allele.tsv")  , emit: tlens
    tuple val(meta), path("${prefix}/*.png")                , emit: plots
    tuple val(meta), path("${prefix}/qc/cmd.txt")           , emit: cmd
    tuple val(meta), path("${prefix}/qc/stats.txt")         , emit: stats
    tuple val(meta), path("${prefix}/qc/qc_readlens.png")   , emit: qc_readlens
    tuple val(meta), path("${prefix}/qc/readlens.npz")      , emit: readlens
    tuple val(meta), path("${prefix}/qc/rng.txt")           , emit: rng
    tuple val("${task.process}"), val('telogator2'), eval("telogator2 --version | sed 's/telogator2 //'"), emit: versions_telogator2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ref_arg = fasta ? "--ref ${fasta}" : ""
    """
    telogator2 \\
        -i ${reads} \\
        -o ${prefix} \\
        -p ${task.cpus} \\
        ${ref_arg} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/qc
    touch ${prefix}/tlens_by_allele.tsv
    touch ${prefix}/all_final_alleles.png
    touch ${prefix}/violin_atl.png
    touch ${prefix}/qc/cmd.txt
    touch ${prefix}/qc/qc_readlens.png
    touch ${prefix}/qc/readlens.npz
    touch ${prefix}/qc/rng.txt
    touch ${prefix}/qc/stats.txt
    """
}
