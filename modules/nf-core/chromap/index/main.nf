process CHROMAP_INDEX {
    tag "$fasta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chromap:0.2.6--hdcf5f25_0' :
        'biocontainers/chromap:0.2.6--hdcf5f25_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path ("*.index"), emit: index
    tuple val("${task.process}"), val('chromap'), eval("chromap --version 2>&1"), topic: versions, emit: versions_chromap

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    chromap \\
        -i \\
        $args \\
        -t $task.cpus \\
        -r $fasta \\
        -o ${prefix}.index
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    touch ${prefix}.index
    """
}
