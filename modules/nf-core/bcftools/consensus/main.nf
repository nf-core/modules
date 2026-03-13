process BCFTOOLS_CONSENSUS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'}"

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
            ${args} \\
            ${masking} \\
            > ${prefix}.fa
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa
    """
}
