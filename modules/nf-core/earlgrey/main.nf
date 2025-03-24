process EARLGREY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/earlgrey:5.1.0--h9948957_0':
        'quay.io/biocontainers/earlgrey:5.1.1--h9948957_0' }"

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${meta.id}"), emit: output_directory
    tuple val(meta), path("${meta.id}/*.fasta"), emit: fasta_files
    tuple val(meta), path("${meta.id}/*.gff"), emit: gff_files
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    earlgrey \\
        $args \\
        -i $genome \\
        -s ${prefix} \\
        -t $task.cpus \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        earlgrey: \$(earlgrey --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    mkdir ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        earlgrey: \$(echo 'v5.1.1')
    END_VERSIONS
    """
}
