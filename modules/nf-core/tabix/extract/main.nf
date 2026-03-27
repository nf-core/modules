process TABIX_EXTRACT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92859404d861ae01afb87e2b789aebc71c0ab546397af890c7df74e4ee22c8dd/data' :
        'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(regions)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"),     emit: vcf, optional: true
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi, optional: true
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'"), topic: versions, emit: versions_tabix

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args  ?: ''
    def args2       = task.ext.args2 ?: ''
    prefix          = task.ext.prefix ?: "${meta.id}"
    def regions_arg = regions ? "-R ${regions}" : ""
    """
    tabix \\
        ${regions_arg} \\
        ${args} \\
        ${vcf} \\
        | bgzip ${args2} > ${prefix}.vcf.gz

    tabix ${prefix}.vcf.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | bgzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
