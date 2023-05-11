process SGDEMUX {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::sgdemux=1.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sgdemux:1.1.1--ha982bd6_0' :
        'biocontainers/sgdemux:1.1.1--ha982bd6_0' }"

    input:
    // Input fastq's must be bgzipped for compatibility with sgdemux
    tuple val(meta), path(sample_sheet), path(fastqs_dir)

    output:
    tuple val(meta), path('output/*_R*.fastq.gz')                   , emit: sample_fastq
    tuple val(meta), path('output/metrics.tsv')                     , emit: metrics
    tuple val(meta), path('output/most_frequent_unmatched.tsv')     , emit: most_frequent_unmatched
    tuple val(meta), path('output/per_project_metrics.tsv')         , emit: per_project_metrics
    tuple val(meta), path('output/per_sample_metrics.tsv')          , emit: per_sample_metrics
    tuple val(meta), path('output/sample_barcode_hop_metrics.tsv')  , emit: sample_barcode_hop_metrics
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p output/
    sgdemux \\
        --sample-metadata ${sample_sheet} \\
        --fastqs ${fastqs_dir} \\
        --output-dir output \\
        --demux-threads ${task.cpus} \\
        --compressor-threads ${task.cpus} \\
        --writer-threads ${task.cpus} \\
        ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sgdemux: \$(echo \$(sgdemux --version 2>&1) | cut -d " " -f2)
    END_VERSIONS
    """
}
