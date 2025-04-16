process COVERM_CONTIG {
    tag "${meta.id}"
    label "process_medium"

    conda "bioconda::coverm=0.7.0-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coverm:0.7.0--h07ea13f_0' :
        'biocontainers/coverm:0.7.0--h07ea13f_0' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(reference)
    val bam_input
    val interleaved

    output:
    tuple val(meta), path("*.depth.txt"), emit: coverage
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ""
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def fastq_input   = meta.single_end ? "--single" : interleaved ? "--interleaved" : "--coupled"
    def input_type    = bam_input ? "--bam-files" : "${fastq_input}"
    def reference_str = bam_input ? "" : "--reference ${reference}"
    """
    TMPDIR=.

    coverm contig \\
        --threads ${task.cpus} \\
        ${input_type} ${input} \\
        ${reference_str} \\
        ${args} \\
        --output-file ${prefix}.depth.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(coverm --version | sed 's/coverm //')
    END_VERSIONS
    """

    stub:
    def prefix        = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.depth.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(coverm --version | sed 's/coverm //')
    END_VERSIONS
    """
}
