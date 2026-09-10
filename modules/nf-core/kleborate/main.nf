process KLEBORATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kleborate:2.1.0--pyhdfd78af_1' :
        'quay.io/biocontainers/kleborate:2.1.0--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('kleborate'), eval("kleborate --version | sed 's/Kleborate v//'"), emit: versions_kleborate, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kleborate \\
        $args \\
        --outfile ${prefix}.results.txt \\
        --assemblies $fastas

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.results.txt

    """
}
