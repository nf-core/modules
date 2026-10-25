process BCFTOOLS_CONSENSUS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/164b5421d87e69f8e9f5add5a937e148e8165f11bfb18159da607f56ecdb1495/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--4e7d12da6416a855'}"

    input:
    tuple val(meta), path(vcf), path(tbi), path(fasta), path(mask)

    output:
    tuple val(meta), path('*.fa'), emit: fasta
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def masking = mask ? "-m ${mask}" : ""
    """
    cat ${fasta} \\
        | bcftools \\
            consensus \\
            ${vcf} \\
            ${masking} \\
            ${args} \\
            > ${prefix}.fa
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa
    """
}
