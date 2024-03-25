process TOULLIGQC {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/toulligqc:2.5.4--pyhdfd78af_0':
        'biocontainers/toulligqc:2.5.4--pyhdfd78af_0' }"

    input:
    path seq_summary
    path fastq
    path bam

    output:
    path "*/*.data"                 , emit: report_data
    path "*/*.html"                 , emit: report_html, optional: true
    path "*/images/*.html"          , emit: plots_html
    path "*/images/plotly.min.js"   , emit: plotly_js

    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def seq_summary_arg = seq_summary ? "--sequencing-summary-source ${seq_summary}" : ""
    def fastq_arg = fastq ? "--fastq ${fastq}" : ""
    def bam_arg = bam ? "--bam ${bam}" : ""

    """
    toulligqc ${seq_summary_arg} \\
                ${fastq_arg} ${bam_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        toulligqc: \$(toulligqc --version 2>&1)
    END_VERSIONS
    """
}
