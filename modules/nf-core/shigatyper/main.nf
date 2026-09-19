process SHIGATYPER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shigatyper:2.0.5--pyhdfd78af_0':
        'quay.io/biocontainers/shigatyper:2.0.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.tsv")     , emit: tsv
    tuple val(meta), path("${prefix}-hits.tsv"), optional: true, emit: hits
    tuple val("${task.process}"), val('shigatyper'), eval("shigatyper --version | sed 's/ShigaTyper //'"), emit: versions_shigatyper, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        """
        shigatyper \\
            ${args} \\
            --SE ${reads} \\
            --name ${prefix}
        """
    } else {
        """
        shigatyper \\
            ${args} \\
            --R1 ${reads[0]} \\
            --R2 ${reads[1]} \\
            --name ${prefix}
        """
    }

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
