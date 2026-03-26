process CUTESV {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/78/78322fdee2a195b18a56432f8b3bbc92b75015b6c921b364d82f0655461992f5/data' :
        'community.wave.seqera.io/library/cutesv:2.1.3--858bd1cbe0a6dc2f' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val("cuteSV"), eval("cuteSV --version 2>&1 | sed 's/cuteSV //'"), topic: versions, emit: versions_cutesv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cuteSV \
        ${bam} \\
        ${fasta} \\
        ${prefix}.vcf \\
        . \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    """
}
