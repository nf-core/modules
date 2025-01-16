process OPENMS_PEPTIDEINDEXER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.2.0--haddbca4_4':
        'biocontainers/openms:3.2.0--haddbca4_4' }"

    input:
    tuple val(meta), path(id_file)
    tuple val(meta), path(id_fasta)

    output:
    tuple val(meta), path("*.idXML"), emit: id_file_pi
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_pi"

    """
    PeptideIndexer \\
        -in $id_file \\
        -fasta $id_fasta \\
        -out ${prefix}.idXML \\
        -threads $task.cpus \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_pi"

    """
    touch ${prefix}.idXML

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """
}
