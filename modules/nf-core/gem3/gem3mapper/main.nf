process GEM3_GEM3MAPPER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-240a9c1936dd6a68f46aa198b2629b6734a18428:543223d3cc2f69d86af72f7e9a3200812ae25327-0':
        'quay.io/biocontainers/mulled-v2-240a9c1936dd6a68f46aa198b2629b6734a18428:543223d3cc2f69d86af72f7e9a3200812ae25327-0' }"

    input:
    tuple val(meta), path(gem)
    tuple val(meta2), path(fastq)
    val   sort_bam

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('gem3-mapper'), eval("gem-mapper --version 2>&1 | sed 's/v//'"), emit: versions_gem3, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    """
    gem-mapper \\
        -F 'SAM' \\
        -I ${gem} \\
        -i ${fastq} \\
        -t ${task.cpus} \\
        ${args} \\
        | samtools ${samtools_command} \\
        ${args2} \\
        -@ ${task.cpus} \\
        -o ${prefix}.bam \\
        -
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
