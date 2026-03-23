process TELOMEREHUNTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/03/032cd65d6d2c4637b637f1f7547d073d097e58fb41a5fb3c6189a94e1df7f217/data'
        : 'community.wave.seqera.io/library/python_pysam_samtools_numpy_pruned:db6e043542c8de22' }"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(control_bam), path(control_bai)

    output:
    tuple val(meta), path("${prefix}/${prefix}_summary.tsv")       , emit: summary
    tuple val(meta), path("${prefix}/tumor_TelomerCnt_${prefix}/")  , emit: tumor
    tuple val(meta), path("${prefix}/control_TelomerCnt_${prefix}/"), emit: control, optional: true
    tuple val("${task.process}"), val('telomerehunter'), eval("pip show telomerehunter 2>/dev/null | grep Version | cut -d ' ' -f2"), emit: versions_telomerehunter, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def control = control_bam ? "-ibc ${control_bam}" : ""
    """
    telomerehunter \\
        -ibt ${tumor_bam} \\
        ${control} \\
        -o . \\
        -p ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/tumor_TelomerCnt_${prefix}
    mkdir -p ${prefix}/control_TelomerCnt_${prefix}
    touch ${prefix}/${prefix}_summary.tsv
    touch ${prefix}/tumor_TelomerCnt_${prefix}/${prefix}_tumor_summary.tsv
    touch ${prefix}/control_TelomerCnt_${prefix}/${prefix}_control_summary.tsv
    """
}
