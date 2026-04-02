process COVERM_CONTIG {
    tag "${meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c4/c4402aed2ff8bd062a013369c9049e68bdd159e849b271a227a8ce2727e941e1/data'
        : 'community.wave.seqera.io/library/coverm:0.7.0--c57cbb9db4d660fb'}"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(reference)
    val bam_input
    val interleaved
    val enable_bam_output

    output:
    tuple val(meta), path('*.depth.tsv'), emit: coverage
    tuple val(meta), path('*.bam')      , emit: bam, optional: true
    tuple val("${task.process}"), val('coverm'), eval('coverm --version | sed "s/coverm //"'), emit: versions_coverm, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args ?: ""
    def prefix         = task.ext.prefix ?: "${meta.id}"
    def fastq_input    = meta.single_end ? "--single" : interleaved ? "--interleaved" : "--coupled"
    def input_type     = bam_input ? "--bam-files" : "${fastq_input}"

    def reference_str  = bam_input ? "" : "--reference ${reference}"
    def bam_output_str = enable_bam_output ? "--bam-file-cache-directory _bam_cache/" : ""
    """
    TMPDIR=.

    coverm contig \\
        --threads ${task.cpus} \\
        ${input_type} ${input} \\
        ${reference_str} \\
        ${bam_output_str} \\
        ${args} \\
        --output-file ${prefix}.depth.tsv

    mv _bam_cache/*.bam . || true
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.depth.tsv
    touch ${prefix}.bam
    """
}
