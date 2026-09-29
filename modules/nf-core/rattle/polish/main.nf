process RATTLE_POLISH {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/rattle:1.0--h5ca1c30_0'
        : 'quay.io/biocontainers/rattle:1.0--h5ca1c30_0'}"

    input:
    tuple val(meta), path(consensi)

    output:
    tuple val(meta), path("${prefix}.transcriptome.fq")  , emit: transcriptome
    tuple val(meta), path("${prefix}.polish_summary.tsv"), emit: summary, optional: true
    tuple val(meta), path("${prefix}.log")               , emit: log
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    tuple val("${task.process}"), val('rattle'), val('1.0'), emit: versions_rattle, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    rattle polish \\
        ${args} \\
        -i ${consensi} \\
        -o . \\
        -t ${task.cpus} \\
        2>| >(tee ${prefix}.log >&2)

    mv transcriptome.fq ${prefix}.transcriptome.fq
    if [[ -f polish_summary.tsv ]]; then
        mv polish_summary.tsv ${prefix}.polish_summary.tsv
    fi
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def summary = args.contains('--summary') ? "printf 'transcript_cluster_0, gene_cluster_0, new_cluster_0\\n' > ${prefix}.polish_summary.tsv" : ''
    """
    printf '@transcript_cluster_0 gene_cluster_0 generated_from_transcript_clusters=1 total_reads=1 labels=\\nACGT\\n+\\nIIII\\n' > ${prefix}.transcriptome.fq
    touch ${prefix}.log
    ${summary}
    """
}
