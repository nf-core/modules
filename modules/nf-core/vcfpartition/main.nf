process VCFPARTITION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/77/7713d869c8c8259c10701c95fc105bad8bcfbd6735de6941a47c9e6e26e9bb2f/data':
        'community.wave.seqera.io/library/bio2zarr:0.1.8--c2c92dd3f64fb0f9' }"

    input:
    tuple val(meta), path(vcf), path(index)

    output:
    tuple val(meta), path("*.tsv"), emit: partitions
    tuple val("${task.process}"), val('vcfpartition'), eval("vcfpartition --version | sed 's/.* //'"), topic: versions, emit: versions_vcfpartition

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vcfpartition \\
        ${args} \\
        ${vcf} \\
        > ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
