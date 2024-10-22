process TIDK_SEARCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tidk:0.2.41--hdbdd923_0':
        'biocontainers/tidk:0.2.41--hdbdd923_0' }"

    input:
    tuple val(meta), path(fasta)
    val string

    output:
    tuple val(meta), path("*.tsv")          , emit: tsv         , optional: true
    tuple val(meta), path("*.bedgraph")     , emit: bedgraph    , optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tidk \\
        search \\
        --string $string \\
        --output $prefix \\
        --dir tidk \\
        $args \\
        $fasta

    mv \\
        tidk/${prefix}_telomeric_repeat_windows.tsv \\
        ${prefix}.tsv \\
        || echo "TSV file was not produced"

    mv \\
        tidk/${prefix}_telomeric_repeat_windows.bedgraph \\
        ${prefix}.bedgraph \\
        || echo "BEDGRAPH file was not produced"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tidk: \$(tidk --version | sed 's/tidk //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--extension bedgraph") ? 'bedgraph' : 'tsv'
    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tidk: \$(tidk --version | sed 's/tidk //')
    END_VERSIONS
    """
}
