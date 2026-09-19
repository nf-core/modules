process GENMAP_INDEX {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genmap:1.3.0--h1b792b2_1' :
        'quay.io/biocontainers/genmap:1.3.0--h1b792b2_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}") , emit: index
    tuple val("${task.process}"), val('genmap'), eval("genmap --version |& sed -n 's/GenMap version: //p'"), emit: versions_genmap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "$meta.id"

    """
    genmap \\
        index \\
        --fasta-file ${fasta} \\
        --index ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "$meta.id"

    """
    touch ${prefix}
    """
}
