process MAKO {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f8/f8581689e70d48ca3b393654c08f87e95a7bed8fb825913dc5d89822c0a4bf4c/data'
:         'community.wave.seqera.io/library/fg-mako:0.1.8--55317497dabc45f3' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam"),     emit: bam
    tuple val(meta), path("${prefix}.bam.bai"), emit: bai, optional: true
    tuple val("${task.process}"), val('mako'),  eval("mako --version | sed '1!d; s/^[^ ]* //; s/ .*//'"),    topic: versions, emit: versions_mako
    tuple val("${task.process}"), val('fgumi'), eval("mako --version | sed '2!d; s/.* //'"),    topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${prefix}.bam" == "${bam}") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    mako \\
        --input ${bam} \\
        --output ${prefix}.bam \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def create_index = args.contains('--write-index') ? "touch ${prefix}.bam.bai" : ''
    """
    touch ${prefix}.bam
    ${create_index}
    """
}
