process GAPSEQ_DOALL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gapseq:1.4.0--hdfd78af_0' :
        'biocontainers/gapseq:1.4.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(medium)

    output:
    tuple val(meta), path("*.RDS")  , emit: model
    tuple val(meta), path("*.csv")  , emit: csv
    tuple val(meta), path("*.fna")  , emit: fna      , optional: true
    tuple val(meta), path("*.log")  , emit: log      , optional: true
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def medium_arg = medium ? "-m $medium" : ''
    """
    gapseq \\
        doall \\
        -p ${task.cpus} \\
        $medium_arg \\
        $args \\
        $fasta

    # Rename output files with prefix
    for file in *.RDS; do
        if [ -f "\$file" ]; then
            mv "\$file" "${prefix}_\$file"
        fi
    done

    for file in *.csv; do
        if [ -f "\$file" ]; then
            mv "\$file" "${prefix}_\$file"
        fi
    done

    for file in *.fna; do
        if [ -f "\$file" ]; then
            mv "\$file" "${prefix}_\$file"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gapseq: \$(gapseq --version 2>&1 | sed 's/gapseq //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_model-filled.RDS
    touch ${prefix}_pathways.csv
    touch ${prefix}_transporters.csv
    """
}
