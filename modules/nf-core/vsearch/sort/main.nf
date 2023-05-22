process VSEARCH_SORT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::vsearch=2.21.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.21.1--h95f258a_0':
        'biocontainers/vsearch:2.21.1--h95f258a_0' }"

    input:
    tuple val(meta), path(fasta)
    val sort_arg

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$fasta" == "${prefix}.fasta") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    vsearch \\
        $sort_arg $fasta \\
        --threads $task.cpus \\
        --output ${prefix}.fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //g;s/,.*//g;s/^v//;s/_.*//')
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //g;s/,.*//g;s/^v//;s/_.*//')
    END_VERSIONS
    """
}
