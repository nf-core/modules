process BEDOPS_CONVERT2BED {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.42--h9948957_0':
        'biocontainers/bedops:2.4.42--h9948957_0' }"

    input:
    tuple val(meta), path(in_file)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val("bedops"), eval('convert2bed --version | sed -n "s/.*version: *\\([^ ]*\\).*/\\1/p"'), emit: versions_bedops, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = in_file.getExtension()
    """
    convert2bed \\
        $args \\
        -i $format \\
        < $in_file \\
        > ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    """
}
