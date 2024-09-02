process MCRONI {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mcroni%3A1.0.4--pyh5e36f6f_0':
        'biocontainers/mcroni:1.0.4--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv")               , emit: tsv
    tuple val(meta), path("*.fa"), optional: true, emit: fa
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    def VERSION = '1.0.4' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    mcroni \\
        $args \\
        --output $prefix \\
        --fasta $fasta_name

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mcroni: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.4' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.tsv
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mcroni: $VERSION
    END_VERSIONS
    """
}
