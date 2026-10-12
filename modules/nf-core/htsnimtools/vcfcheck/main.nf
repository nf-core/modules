process HTSNIMTOOLS_VCFCHECK {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hts-nim-tools:0.3.11--hbeb723e_0':
        'quay.io/biocontainers/hts-nim-tools:0.3.11--hbeb723e_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(background_vcf), path(background_tbi)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('htsnimtools'), eval("hts_nim_tools | sed -n 's/version: //p'"), emit: versions_htsnimtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hts_nim_tools \\
        vcf-check \\
        ${args} \\
        ${background_vcf} \\
        ${vcf} \\
        > ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
