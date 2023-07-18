process MASH_SKETCH {
    tag "$meta.id"
    label 'process_medium'
    conda "bioconda::mash=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.3--he348c14_1' :
        'biocontainers/mash:2.3--he348c14_1' }"

    input:
    tuple val(meta), path(input_seq)
    val(read_sketch)
    val(multi_fasta_sketch)

    output:
    tuple val(meta), path("*.msh")        , emit: mash
    tuple val(meta), path("*.mash_stats") , emit: stats
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_arg = read_sketch ? '-r' : ''
    def multi_fasta_arg = multi_fasta_sketch ? '-i' : ''
    """
    mash \\
        sketch \\
        $args \\
        -p $task.cpus \\
        -o ${prefix} \\
        $read_arg \\
        $multi_fasta_arg \\
        $input_seq \\
        2> ${prefix}.mash_stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_arg = read_sketch ? '-r' : ''
    def multi_fasta_arg = multi_fasta_sketch ? '-i' : ''
    """
    touch ${prefix}.msh
    touch ${prefix}.mash_stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1)
    END_VERSIONS
    """
}
