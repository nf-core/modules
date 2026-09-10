process TD2_LONGORFS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ec/ecc5fb9460ed4255919340ea296994a3491c91515e9583c842ed49599e22a9ce/data':
        "community.wave.seqera.io/library/td2:1.1.0--76046a413a4219c1"}"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), val(min_len), val(abs_min_len)
    tuple val(meta3), path(gene_transmap)

    output:
    tuple val(meta), path("${prefix}/longorfs/longest_orfs.{cds,pep,gff3}"), emit: orfs
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('TD2.LongOrfs'), val("1.1.0"), emit: versions_td2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def min_length = min_len ?: '90'
    def absolute_min_length = abs_min_len ?: '90'
    def gene_trans_map = gene_transmap ? "--gene-trans-map ${gene_transmap}" : ''

    """
    TD2.LongOrfs \\
        -t ${fasta} \\
        --min-length ${min_length} \\
        --absolute-min-length ${absolute_min_length} \\
        ${gene_trans_map} \\
        --output-dir ${prefix}/longorfs \\
        --threads ${task.cpus} \\
        --memory-threshold ${task.memory.toGiga()} \\
        --verbose \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}/longorfs
    touch ${prefix}/longorfs/longest_orfs.cds
    touch ${prefix}/longorfs/longest_orfs.gff3
    touch ${prefix}/longorfs/longest_orfs.pep
    """
}
