process BCFTOOLS_PLUGINVCF2TABLE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9e674cbc269b944e50ad4f5846ee92a2196b59d9ed76e9146867c76fb32cb8aa/data'
        : 'community.wave.seqera.io/library/bcftools:1.23.1--4d193a5f61d4aed7'}"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(regions), path(targets), path(samples)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions_file = regions ? "--regions-file ${regions}" : ''
    def targets_file = targets ? "--targets-file ${targets}" : ''
    def samples_file = samples ? "--samples-file ${samples}" : ''
    """
    bcftools \\
        +vcf2table \\
        ${args} \\
        ${regions_file} \\
        ${targets_file} \\
        ${samples_file} \\
        ${vcf} \\
        > ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
