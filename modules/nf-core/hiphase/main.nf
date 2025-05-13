process HIPHASE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hiphase:1.4.5--h9ee0642_0':
        'biocontainers/hiphase:1.4.5--h9ee0642_0' }"

    input:
    tuple val(meta), path(vcf), path(csi)
    tuple val(meta2), path(bam), path(bai)
    tuple val(meta3), path(fasta)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(meta), path("*.csv"), emit: csv
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
        --output-vcf ${prefix}.phased.vcf \
        --output-bam ${prefix}.phased.bam \
        --reference $fasta \
        --stats-file ${prefix}.stats.csv \
        $args \
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hiphase: \$(hiphase --version |& sed '1!d ; s/hiphase //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.phased.vcf
    touch ${prefix}.stats.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hiphase: \$(hiphase --version |& sed '1!d ; s/hiphase //')
    END_VERSIONS
    """
}
