process SVIM_ALIGNMENT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svim:2.0.0--pyhdfd78af_0':
        'biocontainers/svim:2.0.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(task.process), val('svim'), val('2.0.0'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Set matplotlib config to avoid warnings
    export MPLCONFIGDIR=\$(mktemp -d)

    svim alignment \\
        ${args} \\
        ${prefix} \\
        ${bam} \\
        ${fasta}

    mv ${prefix}/variants.vcf ${prefix}.vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.vcf"
    """
}