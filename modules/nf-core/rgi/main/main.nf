process RGI_MAIN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::rgi=6.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rgi:6.0.2--pyha8f3691_0':
        'biocontainers/rgi:6.0.2--pyha8f3691_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.txt") , emit: tsv
    tuple val(meta), path("temp/")  , emit: tmp
    env VER                        , emit: tool_version
    env DBVER                      , emit: db_version
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    rgi \\
        main \\
        $args \\
        --num_threads $task.cpus \\
        --output_file $prefix \\
        --input_sequence $fasta

    mkdir temp/
    mv *.xml *.fsa *.{nhr,nin,nsq} *.draft *.potentialGenes *{variant,rrna,protein,predictedGenes,overexpression,homolog}.json temp/

    VER=\$(rgi main --version)
    DBVER=\$(rgi database --version)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(echo \$VER)
        rgi-database: \$(echo \$DBVER)
    END_VERSIONS
    """
}
