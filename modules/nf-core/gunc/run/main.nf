process GUNC_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gunc:1.0.6--pyhdfd78af_0' :
        'biocontainers/gunc:1.0.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta_files, stageAs: 'input_files/*')
    path(db)

    output:
    tuple val(meta), path("*maxCSS_level.tsv")                  , emit: maxcss_level_tsv
    tuple val(meta), path("*all_levels.tsv")    , optional: true, emit: all_levels_tsv
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ls input_files/* > input_files.txt
    gunc \\
        run \\
        --input_file input_files.txt \\
        --db_file $db \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunc: \$( gunc --version )
    END_VERSIONS
    """

    stub:
    """
    touch maxCSS_level.tsv all_levels.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunc: \$( gunc --version )
    END_VERSIONS
    """
}
