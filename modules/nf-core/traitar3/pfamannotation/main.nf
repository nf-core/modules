process TRAITAR3_PFAMANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hdbdd923_1' :
        'biocontainers/hmmer:3.4--hdbdd923_1' }"

    input:
    tuple val(meta), path(faa)
    path pfam_hmm

    output:
    tuple val(meta), path("*.domtblout.gz"), emit: domtblout
    tuple val(meta), path("*.txt.gz")      , emit: output
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hmmsearch \\
        --notextw \\
        --cpu ${task.cpus} \\
        --cut_ga \\
        --domtblout ${prefix}.domtblout \\
        ${args} \\
        -o ${prefix}.txt \\
        ${pfam_hmm} \\
        ${faa}

    gzip --no-name ${prefix}.domtblout ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmsearch -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.domtblout
    touch ${prefix}.txt

    gzip --no-name ${prefix}.domtblout ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmsearch -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """
}
