process BIOBAMBAM_BAMMARKDUPLICATES2 {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/biobambam:2.0.185--h85de650_1'
        : 'biocontainers/biobambam:2.0.185--h85de650_1'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.metrics.txt"), emit: metrics
    tuple val("${task.process}"), val('biobambam'), eval("bammarkduplicates2 --version |& sed '1!d; s/.*version //; s/.\$//'"), topic: versions, emit: versions_biobambam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bammarkduplicates2 \\
        ${args} \\
        I=${bam} \\
        O=${prefix}.bam \\
        M=${prefix}.metrics.txt \\
        tmpfile=${prefix} \\
        markthreads=${task.cpus}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.metrics.txt
    """
}
