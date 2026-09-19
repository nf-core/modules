process HAPLOGREP2_CLASSIFY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/haplogrep:2.4.0--hdfd78af_0':
        'quay.io/biocontainers/haplogrep:2.4.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(inputfile)
    val(format)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('haplogrep2'), eval("haplogrep --version 2>&1 | sed '2!d;s/.*v//'"), emit: versions_haplogrep2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    haplogrep \\
        classify \\
        ${args} \\
        --in ${inputfile} \\
        --out ${prefix}.txt \\
        --format ${format}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """

}
