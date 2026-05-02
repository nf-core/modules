process SLIMFASTQ {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/slimfastq:2.04--h87f3376_2':
        'quay.io/biocontainers/slimfastq:2.04--h87f3376_2' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.sfq"), emit: sfq
    tuple val("${task.process}"), val('slimfastq'), eval('echo 2.04'), emit: versions_slimfastq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        gzip -d -c '${fastq}' | slimfastq \\
            $args \\
            -f '${prefix}.sfq'
        """
    } else {
        """
        gzip -d -c '${fastq[0]}' | slimfastq \\
            $args \\
            -f '${prefix}_1.sfq'

        gzip -d -c '${fastq[1]}' | slimfastq \\
            $args \\
            -f '${prefix}_2.sfq'
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sfq
    """
}
