process TOULLIGQC {
    label 'process_low'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/toulligqc:2.5.4--pyhdfd78af_0':
        'biocontainers/toulligqc:2.5.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(seq_summary)
    tuple val(meta), path(fastq)
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*/*.data")   , emit: report_data
    path "*/*.html"                     , emit: report_html, optional: true
    path "*/images/*.html"              , emit: plots_html
    path "*/images/plotly.min.js"       , emit: plotly_js

    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def seq_summary_input = seq_summary ? "--sequencing-summary-source ${seq_summary}" : ""
    def fastq_input = fastq ? "--fastq ${fastq}" : ""
    def bam_input = bam ? "--bam ${bam}" : ""

    """
    toulligqc ${seq_summary_input} \\
                ${fastq_input} \\
                ${bam_input} \\
                $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        toulligqc: \$(toulligqc --version 2>&1)
    END_VERSIONS
    """
}
