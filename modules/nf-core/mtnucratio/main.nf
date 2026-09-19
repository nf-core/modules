process MTNUCRATIO {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mtnucratio:0.7--hdfd78af_2' :
        'quay.io/biocontainers/mtnucratio:0.7--hdfd78af_2' }"

    input:
    tuple val(meta), path(bam)
    val(mt_id)

    output:
    tuple val(meta), path("*.mtnucratio"), emit: mtnucratio
    tuple val(meta), path("*.json")      , emit: json
    tuple val("${task.process}"), val('mtnucratio'), eval("mtnucratio --version 2>&1 | sed -n 's/Version: //p'"), emit: versions_mtnucratio, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mtnucratio \\
        ${args} \\
        ${bam} \\
        ${mt_id}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mtnucratio
    touch ${prefix}.json
    """
}
