process YARA_INDEX {
    tag "${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::yara=1.0.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/yara:1.0.2--2' :
        'quay.io/biocontainers/yara:1.0.2--2' }"

    input:
    val(meta)
    path fasta

    output:
    tuple val(meta), path("${prefix}*")  , emit: index
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"

    """
    yara_indexer \\
        $fasta \\
        -o $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yara: \$(echo \$(yara_indexer --version 2>&1) | sed 's/^.*yara_indexer version: //; s/ .*\$//')
    END_VERSIONS
    """
}
