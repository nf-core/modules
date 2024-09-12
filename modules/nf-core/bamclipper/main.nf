
process BAMCLIPPER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamclipper:1.0.0--hdfd78af_2':
        'biocontainers/bamclipper:1.0.0--hdfd78af_2' }"

    input:
    tuple val(meta), path(bam), path(bai), path(bedpe)

    output:
    tuple val(meta), path("*.primerclipped.bam")    , emit: bam
    tuple val(meta), path("*.primerclipped.bam.bai"), emit: bai
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.0' // WARN: Version information not provided by tool on CLI
    """
    bamclipper.sh \\
        -b ${bam} \\
        -p ${bedpe} \\
        -n $task.cpus \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamclipper: $VERSION
        samtools: \$( samtools --version |& sed '1!d ; s/samtools //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.0'
    """
    touch ${prefix}.primerclipped.bam
    touch ${prefix}.primerclipped.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamclipper: $VERSION
        samtools: \$( samtools --version |& sed '1!d ; s/samtools //' )
    END_VERSIONS
    """
}
