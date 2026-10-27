process BCFTOOLS_PLUGINVCF2TABLE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/97/97495e3fd16a7177cfaa1fa318556ce3bbe372d1234c94ee5b9f78d7fde6777c/data'
        : 'community.wave.seqera.io/library/bcftools:1.23.1--e6b944fbe2dde7a7'}"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(regions), path(targets), path(samples)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

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
