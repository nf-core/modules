
process BAMCLIPPER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamclipper:1.0.0--hdfd78af_2':
        'biocontainers/bamclipper:1.0.0--hdfd78af_2' }"

    input:
    tuple val(meta), path(bam), path(bai), path(bedpe)

    output:
    tuple val(meta), path("*.primerclipped.bam")    , emit: bam
    tuple val(meta), path("*.primerclipped.bam.bai"), emit: bai
    tuple val("${task.process}"), val('bamclipper'), val("1.0.0"), emit: versions_bamclipper, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions
    // WARN: Version information not provided by tool on CLI
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bamclipper.sh \\
        -b ${bam} \\
        -p ${bedpe} \\
        -n $task.cpus \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.primerclipped.bam
    touch ${prefix}.primerclipped.bam.bai
    """
}
