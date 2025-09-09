process COVERM_GENOME {
    tag "${meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coverm:0.7.0--hcb7b614_4' :
        'biocontainers/coverm:0.7.0--hcb7b614_4' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(reference)
    val bam_input
    val interleaved
    val ref_mode   // "dir" | "file" | "auto"

    output:
    tuple val(meta), path("*.tsv"), emit: coverage
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def _ref_mode = ref_mode ?: 'auto'
    def _allowed = ['dir','file','auto']
    assert _allowed.contains(_ref_mode) : "Invalid ref_mode='${_ref_mode}'. Allowed: ${_allowed.join(', ')}"

    def args          = task.ext.args ?: ""
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def fastq_input   = meta.single_end ? "--single" : interleaved ? "--interleaved" : "--coupled"
    def input_type    = bam_input ? "--bam-files" : "${fastq_input}"

    def reference_str  = (
        _ref_mode == 'dir'  || (_ref_mode == 'auto' && reference.isDirectory())
        ) ? "--genome-fasta-directory ${reference}"
        :   "--genome-fasta-files ${reference}"

    """
    TMPDIR=.

    coverm genome \\
        --threads ${task.cpus} \\
        ${input_type} ${input} \\
        ${reference_str} \\
        ${args} \\
        --output-file ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(coverm --version | sed 's/coverm //')
    END_VERSIONS
    """

    stub:
    def prefix        = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(coverm --version | sed 's/coverm //')
    END_VERSIONS
    """
}
