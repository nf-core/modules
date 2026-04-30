process TIDK_SEARCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tidk:0.2.7--h6872113_0':
        'quay.io/biocontainers/tidk:0.2.7--h6872113_0' }"

    input:
    tuple val(meta), path(fasta)
    val string

    output:
    tuple val(meta), path("*.tsv")          , emit: tsv         , optional: true
    tuple val(meta), path("*.bedgraph")     , emit: bedgraph    , optional: true
    tuple val("${task.process}"), val('tidk'), eval("tidk --version | sed 's/tidk //'"), emit: versions_tidk, topic: versions

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

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--extension bedgraph") ? 'bedgraph' : 'tsv'
    """
    touch ${prefix}.${extension}

    """
}
