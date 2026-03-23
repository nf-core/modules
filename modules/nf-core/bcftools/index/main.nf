process BCFTOOLS_INDEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.csi"), emit: csi, optional: true
    tuple val(meta), path("*.tbi"), emit: tbi, optional: true
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bcftools \\
        index \\
        ${args} \\
        --threads ${task.cpus} \\
        ${vcf}
    """

    stub:
    def args = task.ext.args ?: ''
    def extension = args.contains("--tbi") || args.contains("-t")
        ? "tbi"
        : "csi"
    """
    touch ${vcf}.${extension}
    """
}
