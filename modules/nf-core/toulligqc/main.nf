process TOULLIGQC {
    label 'process_low'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/958e7692f05dbbae6937b082b0e2863f1ec26be7dd6d2b016ae130eaef6cdc85/data'
:         'community.wave.seqera.io/library/toulligqc_huggingface_hub_matplotlib_numpy_pruned:30ed42003347882d' }"

    input:

    tuple val(meta), path(ontfile)

    output:
    tuple val(meta), path("*/*.data"), emit: report_data
    tuple val(meta), path("*/*.html"), emit: report_html, optional: true
    tuple val(meta), path("*/images/*.html"), emit: plots_html
    tuple val(meta), path("*/images/plotly.min.js"), emit: plotly_js
    tuple val("${task.process}"), val('toulligqc'), eval('toulligqc --version'), emit: versions_toulligqc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def input_file = "${ontfile}".endsWith(".fastq") || "${ontfile}".endsWith(".fastq.gz") || "${ontfile}".endsWith(".fq") || "${ontfile}".endsWith(".fq.gz")
        ? "--fastq ${ontfile}"
        : "${ontfile}".endsWith(".txt") || "${ontfile}".endsWith(".txt.gz")
            ? "--sequencing-summary-source ${ontfile}"
            : "${ontfile}".endsWith(".bam") ? "--bam ${ontfile}" : ''

    """
    toulligqc \\
        ${input_file} \\
        --output-directory ${prefix} \\
        ${args}
    """

    stub:
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
    """
}
