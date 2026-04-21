process BCFTOOLS_PLUGINVCF2TABLE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/bcftools:1.23.1--16b1a31e5dc795f7'
        : 'community.wave.seqera.io/library/bcftools:1.23.1--4d193a5f61d4aed7'}"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(regions)
    tuple val(meta3), path(targets)
    tuple val(meta4), path(samples)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

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
