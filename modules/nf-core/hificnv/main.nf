process HIFICNV {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b5/b52ca8757ee3bc9f475e50fa148133ba5600fe1523c1f698347ee63adcb5bba7/data':
        'community.wave.seqera.io/library/hificnv:1.0.1--b7e433ac6789e2d2' }"

    input:
    tuple val(meta) , path(bam), path(bai), path(maf)
    tuple val(meta2), path(ref)
    tuple val(meta3), path(exclude)
    tuple val(meta4), path(expected_cn)

    output:
    tuple val(meta), path("*.copynum.bedgraph"), emit: copynum
    tuple val(meta), path("*.depth.bw")        , emit: depth
    tuple val(meta), path("*.maf.bw")          , emit: maf     , optional: true
    tuple val(meta), path("*.vcf.gz")          , emit: vcf
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // handle optional inputs
    def maf_arg         = maf         ? "--maf ${maf}"                 : ""
    def exclude_arg     = exclude     ? "--exclude ${exclude}"         : ""
    def expected_cn_arg = expected_cn ? "--expected-cn ${expected_cn}" : ""

    """
    hificnv \\
        --bam ${bam} \\
        --ref ${ref} \\
        ${maf_arg} \\
        ${exclude_arg} \\
        ${expected_cn_arg} \\
        --threads ${task.cpus} \\
        --output-prefix ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hificnv: \$(hificnv --version 2>&1 | sed 's/^.*hificnv //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def create_maf = maf ? "touch ${prefix}.maf.bw" : ""

    """
    touch ${prefix}.depth.bw
    touch ${prefix}.copynum.bedgraph
    ${create_maf}
    echo "" | gzip > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hificnv: \$(hificnv --version 2>&1 | sed 's/.* //')
    END_VERSIONS
    """
}
