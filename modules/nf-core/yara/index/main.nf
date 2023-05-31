process YARA_INDEX {
    tag "$fasta"
    label 'process_medium'

    conda "bioconda::yara=1.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/yara:1.0.2--2' :
        'biocontainers/yara:1.0.2--2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta}*")   , emit: index
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''

    """
    yara_indexer \\
        $fasta \\
        -o ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yara: \$(echo \$(yara_indexer --version 2>&1) | sed 's/^.*yara_indexer version: //; s/ .*\$//')
    END_VERSIONS
    """
}
