process SGDEMUX {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sgdemux:1.1.1--ha982bd6_0' :
        'biocontainers/sgdemux:1.1.1--ha982bd6_0' }"

    input:
    // Input fastq's must be bgzipped for compatibility with sgdemux
    tuple val(meta), path(sample_sheet), path(fastqs_dir)

    output:
    tuple val(meta), path("${prefix}/*_R*.fastq.gz")                   , emit: sample_fastq
    tuple val(meta), path("${prefix}/metrics.tsv")                     , emit: metrics
    tuple val(meta), path("${prefix}/most_frequent_unmatched.tsv")     , emit: most_frequent_unmatched
    tuple val(meta), path("${prefix}/per_project_metrics.tsv")         , emit: per_project_metrics
    tuple val(meta), path("${prefix}/per_sample_metrics.tsv")          , emit: per_sample_metrics
    tuple val(meta), path("${prefix}/sample_barcode_hop_metrics.tsv")  , emit: sample_barcode_hop_metrics
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    sgdemux \\
        --sample-metadata ${sample_sheet} \\
        --fastqs ${fastqs_dir} \\
        --output-dir ${prefix} \\
        --demux-threads ${task.cpus} \\
        --compressor-threads ${task.cpus} \\
        --writer-threads ${task.cpus} \\
        ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sgdemux: \$(echo \$(sgdemux --version 2>&1) | cut -d " " -f2)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}

    for sample in \$(cut -d, -f2 ${sample_sheet} | tail -n +2 ); do
        uppercase_sample=\$(echo \$sample | tr '[:lower:]' '[:upper:]')
        echo | gzip > ${prefix}/\${sample}_\${uppercase_sample}_L001_R1_001.fastq.gz
        echo | gzip > ${prefix}/\${sample}_\${uppercase_sample}_L001_R2_001.fastq.gz
    done

    echo | gzip > ${prefix}/out_L001_R1_001.fastq.gz
    echo | gzip > ${prefix}/out_L001_R2_001.fastq.gz

    touch ${prefix}/metrics.tsv
    touch ${prefix}/most_frequent_unmatched.tsv
    touch ${prefix}/per_project_metrics.tsv
    touch ${prefix}/per_sample_metrics.tsv
    touch ${prefix}/sample_barcode_hop_metrics.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sgdemux: \$(echo \$(sgdemux --version 2>&1) | cut -d " " -f2)
    END_VERSIONS
    """
}
