process PBTK_BAM2FASTQ {
    tag "$meta.id"
    label 'process_medium'
    conda "bioconda::pbtk=3.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbtk:3.1.1--h9ee0642_0':
        'biocontainers/pbtk:3.1.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam), path(index)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bam2fastq \\
        $args \\
        -o ${prefix} \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bam2fastx: \$(echo \$(bam2fastq --version 2>&1) | sed 's/^.*bam2fastq //' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bam2fastx: \$(echo \$(bam2fastq --version 2>&1) | sed 's/^.*bam2fastq //' ))
    END_VERSIONS
    """
}
