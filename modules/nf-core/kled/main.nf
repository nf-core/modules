process KLED {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kled:1.2.11--he5eddf3_0' :
        'quay.io/biocontainers/kled:1.2.11--he5eddf3_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val("kled"), eval("kled --version 2>&1 | sed 's/Kled version //' | sed 's/\\.\$//'"), topic: versions, emit: versions_kled

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kled \
        ${args} \\
        --Ref=${fasta} \\
        --threads=${task.cpus} \\
        ${bam} \\
        > ${prefix}.vcf
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    """
}
