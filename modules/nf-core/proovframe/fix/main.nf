process PROOVFRAME_FIX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/proovframe:0.9.7--hdfd78af_1':
        'quay.io/biocontainers/proovframe:0.9.7--hdfd78af_1' }"

    input:
    tuple val(meta) , path(fa)
    tuple val(meta2), path(tsv)

    output:
    tuple val(meta), path("*.fa"), emit: out_fa
    tuple val("${task.process}"), val('proovframe'), eval("proovframe 2>&1 | sed '/proovframe-v/!d;s/proovframe-v//;s/ .*//'"), topic: versions, emit: versions_proovframe

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    proovframe   \\
        fix \\
        ${args} \\
        -o ${prefix}.fa  \\
        ${fa} \\
        ${tsv}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa
    """
}
