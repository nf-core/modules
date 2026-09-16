process CPC2 {
    tag "$meta.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cpc2:1.0.1--hdfd78af_1':
        'quay.io/biocontainers/cpc2:1.0.1--hdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('CPC2'), eval("CPC2.py --version"), topic: versions, emit: versions_cpc2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    CPC2.py \\
       -i ${fasta} \\
       -o ${prefix} \\
    $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    touch ${prefix}.txt
    """
}
