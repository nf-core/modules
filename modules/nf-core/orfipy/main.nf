process ORFIPY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orfipy:0.0.4--py310h184ae93_4':
        'biocontainers/orfipy:0.0.4--py310h184ae93_4' }"

    input:
    tuple val(meta), path(infile)

    output:
    tuple val(meta), path("${prefix}/${prefix}.bed"), emit: bed
    tuple val("${task.process}"), val('orfipy'), eval("orfipy --version | sed 's/.*version //'"), emit: versions_orfipy, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    orfipy \\
        ${infile} \\
        --outdir ${prefix} \\
        --bed ${prefix}.bed \\
        --procs ${task.cpus} \\
        $args
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}/
    touch ${prefix}/${prefix}.bed
    """
}
