process MAPAD_INDEX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mapad:0.45.0--ha96b9cd_0':
        'quay.io/biocontainers/mapad:0.45.0--ha96b9cd_0' }"

    input:
    tuple val(meta), path(fasta, stageAs: "mapad/*")

    output:
    tuple val(meta), path("mapad/"), emit: index
    tuple val("${task.process}"), val("mapad"), eval("mapad --version | sed 's/^mapAD //'"), topic: versions, emit: versions_mapad

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mapad \\
        index \\
        ${args} \\
        --reference ${fasta} \\
        --threads ${task.cpus}
    """

    stub:
    """
    touch ${fasta}.{tbw,tle,toc,tos,tpi,trt,tsa}
    """
}
