process HIPHASE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hiphase:1.4.5--h9ee0642_0':
        'biocontainers/hiphase:1.4.5--h9ee0642_0' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta1), path(index)
    tuple val(meta2), path(bam)
    tuple val(meta3), path(bai)
    tuple val(meta4), path(fasta)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hiphase \
        --bam $bam \
        --vcf $vcf \
        --output-vcf ${prefix}.vcf \
        --output-bam ${prefix}.bam \
        --reference $fasta
        --threads ${task.cpus}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hiphase: \$(hiphase --version |& sed '1!d ; s/hiphase //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hiphase: \$(hiphase --version |& sed '1!d ; s/hiphase //')
    END_VERSIONS
    """
}
