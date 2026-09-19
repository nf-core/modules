process MASHMAP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mashmap:3.1.3--h07ea13f_0':
        'quay.io/biocontainers/mashmap:3.1.3--h07ea13f_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.paf"), emit: paf
    tuple val("${task.process}"), val("mashmap"), eval("mashmap --version 2>&1"), emit: versions_mashmap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mashmap -q ${fasta} \\
        -r ${reference} \\
        -t ${task.cpus} \\
        -o ${prefix}.paf \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.paf
    """
}
