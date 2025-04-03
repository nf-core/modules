process TOULLIGQC {
    label 'process_low'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/toulligqc:2.5.6--pyhdfd78af_0':
        'biocontainers/toulligqc:2.5.6--pyhdfd78af_0' }"

    input:

    tuple val(meta), path(ontfile)


    output:
    tuple val(meta), path("*/*.data")                   , emit: report_data
    tuple val(meta), path("*/*.html")                   , emit: report_html, optional: true
    tuple val(meta), path("*/images/*.html")            , emit: plots_html
    tuple val(meta), path("*/images/plotly.min.js")     , emit: plotly_js
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def input_file = ("$ontfile".endsWith(".fastq") || "$ontfile".endsWith(".fastq.gz") || "$ontfile".endsWith(".fq") || "$ontfile".endsWith(".fq.gz")) ? "--fastq ${ontfile}" :
        ("$ontfile".endsWith(".txt") || "$ontfile".endsWith(".txt.gz")) ? "--sequencing-summary-source ${ontfile}" :
        ("$ontfile".endsWith(".bam")) ? "--bam ${ontfile}" : ''

    """
    toulligqc \\
        $input_file \\
        --output-directory ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        toulligqc: \$(toulligqc --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    mkdir ${prefix}/images
    touch ${prefix}/report.data
    touch ${prefix}/images/Correlation_between_read_length_and_PHRED_score.html
    touch ${prefix}/images/Distribution_of_read_lengths.html
    touch ${prefix}/images/PHRED_score_density_distribution.html
    touch ${prefix}/images/Read_count_histogram.html
    touch ${prefix}/images/plotly.min.js

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        toulligqc: \$(toulligqc --version)
    END_VERSIONS
    """
}
