process BIOBAMBAM_BAMMERGE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/biobambam:2.0.185--h85de650_1'
        : 'biocontainers/biobambam:2.0.185--h85de650_1'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val(meta), path("*.bai"), optional: true, emit: bam_index
    tuple val(meta), path("*.md5"), optional: true, emit: checksum
    tuple val("${task.process}"), val('biobambam'), eval("bammerge --version |& sed '1!d; s/.*version //; s/.\$//'"), topic: versions, emit: versions_biobambam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input_string = bam.join(" I=")

    """
    bammerge \\
        I=${input_string} \\
        ${args} \\
        > ${prefix}.bam
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bam
    """
}
