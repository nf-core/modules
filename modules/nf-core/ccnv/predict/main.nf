process CCNV_PREDICT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ccnv:1.0.0':
        'quay.io/biocontainers/ccnv:1.0.0' }"

    input:
    tuple val(meta), path(reads), path(reads_index), path(ccnv)
    tuple val(meta2), path(ccnv_reference)
    tuple val(meta3), path(fasta), path(fai)
    tuple val(meta4), path(mappability)
    tuple val(meta5), path(gc_content), path(gc_content_index)

    output:
    tuple val(meta), path("*_qc.json"), emit: qc_json
    tuple val(meta), path("*_plot.html"), emit: plot_html, optional: true
    tuple val(meta), path("*_bins.bed.gz"), path("*_bins.bed.gz.tbi"), emit: bins_bed
    tuple val(meta), path("*_segments.bed.gz"), path("*_segments.bed.gz.tbi"), emit: segments_bed
    tuple val(meta), path("*_calls.bed.gz"), path("*_calls.bed.gz.tbi"), emit: calls_bed
    tuple val(meta), path("*_bins.vcf.gz"), path("*_bins.vcf.gz.tbi"), emit: bins_vcf, optional: true
    tuple val(meta), path("*_segments.vcf.gz"), path("*_segments.vcf.gz.tbi"), emit: segments_vcf, optional: true
    tuple val(meta), path("*_calls.vcf.gz"), path("*_calls.vcf.gz.tbi"), emit: calls_vcf, optional: true
    tuple val("${task.process}"), val('ccnv'), eval("ccnv --version 2>&1 | sed 's/ccnv //'"), topic: versions, emit: versions_ccnv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_cmd = ccnv ?: reads
    def fasta_cmd = fasta ? "--fasta $fasta" : ""
    def mappability_cmd = mappability ? "--mappability $mappability" : ""
    def gc_cmd = gc_content ? "--gc $gc_content" : ""

    """
    ccnv \\
        $args \\
        --prefix $prefix \\
        --threads $task.cpus \\
        $gc_cmd \\
        $mappability_cmd \\
        $fasta_cmd \\
        $input_cmd \\
        $ccnv_reference
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}_qc.json
    echo "" | gzip > ${prefix}_bins.bed.gz
    touch ${prefix}_bins.bed.gz.tbi
    echo "" | gzip >  ${prefix}_calls.bed.gz
    touch ${prefix}_calls.bed.gz.tbi
    echo "" | gzip > ${prefix}_calls.vcf.gz
    touch ${prefix}_calls.vcf.gz.tbi
    touch ${prefix}_plot.html
    touch ${prefix}_qc.json
    echo "" | gzip >  ${prefix}_segments.bed.gz
    touch ${prefix}_segments.bed.gz.tbi
    echo "" | gzip > ${prefix}_segments.vcf.gz
    touch ${prefix}_segments.vcf.gz.tbi
    """
}
