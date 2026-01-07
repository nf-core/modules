process GAPSEQ_FIND {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gapseq:1.4.0--hdfd78af_0' :
        'biocontainers/gapseq:1.4.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.csv")  , emit: csv
    tuple val(meta), path("*.fna")  , emit: fna      , optional: true
    tuple val(meta), path("*.log")  , emit: log      , optional: true
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gapseq \\
        find \\
        -p all \\
        -b 200 \\
        -m Bacteria \\
        -t ${task.cpus} \\
        $args \\
        $fasta

    # Rename output files with prefix
    for file in *.csv; do
        if [ -f "\$file" ]; then
            mv "\$file" "${prefix}_\$file"
        fi
    done

    # Optionally rename fna files if they exist
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
    touch ${prefix}_pathways.csv
    touch ${prefix}_reactions.csv
    """
}
