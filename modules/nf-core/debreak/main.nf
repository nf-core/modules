process DEBREAK {
    tag "$meta.id"
    label 'process_high'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/debreak:1.3--h9ee0642_0':
        'quay.io/biocontainers/debreak:1.3--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}/*.vcf"), emit: vcf
    tuple val("${task.process}"), val('debreak'), val('1.3'), emit: versions_debreak, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args  = task.ext.args ?: ''
    def ref_arg = fasta ? "-r ${fasta}" : ""
    """
    debreak \\
        --bam ${bam} \\
        -o ${prefix} \\
        -t ${task.cpus} \\
        ${ref_arg} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/debreak.vcf
    """
}
