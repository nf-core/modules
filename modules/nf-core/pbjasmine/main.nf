process PBJASMINE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbjasmine:26.1.3--hd63eeec_0':
        'quay.io/biocontainers/pbjasmine:26.1.3--hd63eeec_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('pbjasmine'), eval("jasmine --version | head -n 1 | sed 's/jasmine //'"), topic: versions, emit: versions_pbjasmine


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_jasmine"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    jasmine \\
        $args \\
        --num-threads $task.cpus \\
        $bam \\
        ${prefix}.bam

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_jasmine"

    """
    touch ${prefix}.bam
    """
}
