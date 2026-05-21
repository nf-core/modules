process CNVPYTOR_IMPORTREADDEPTH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvpytor:1.3.1--pyhdfd78af_1':
        'quay.io/biocontainers/cnvpytor:1.3.1--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(input_file), path(index)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.pytor"), emit: pytor
    tuple val("${task.process}"), val('cnvpytor'), eval("cnvpytor --version 2>&1 | sed -n 's/.*CNVpytor //p'"), emit: versions_cnvpytor, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "-T ${fasta}" : ''
    """
    cnvpytor \\
        -root ${prefix}.pytor \\
        -rd $input_file \\
        $args \\
        $reference
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pytor
    """
}
