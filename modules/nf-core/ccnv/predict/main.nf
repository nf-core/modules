process CCNV_PREDICT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cd/cd552e195217d37801cf86348dbe875fa66ae2014147e481fd6be6e80dbbcd90/data'
        : 'community.wave.seqera.io/library/ccnv:1.1.0--726d8b142feb0e3d'}"

    input:
    tuple val(meta), path(reads), path(reads_index), path(ccnv)
    tuple val(meta2), path(ccnv_reference)
    tuple val(meta3), path(fasta), path(fai)
    tuple val(meta4), path(mappability)
    tuple val(meta5), path(gc_content), path(gc_content_index)

    output:
    tuple val(meta), path("*_qc.json"), emit: qc_json
    tuple val(meta), path("*_plot.html"), emit: plot_html, optional: true
    tuple val(meta), path("*_bins.bed.gz"), path("*_bins.bed.gz.{tbi,csi}"), emit: bins_bed
    tuple val(meta), path("*_segments.bed.gz"), path("*_segments.bed.gz.{tbi,csi}"), emit: segments_bed
    tuple val(meta), path("*_calls.bed.gz"), path("*_calls.bed.gz.{tbi,csi}"), emit: calls_bed
    tuple val(meta), path("*_bins.vcf.gz"), path("*_bins.vcf.gz.{tbi,csi}"), emit: bins_vcf, optional: true
    tuple val(meta), path("*_segments.vcf.gz"), path("*_segments.vcf.gz.{tbi,csi}"), emit: segments_vcf, optional: true
    tuple val(meta), path("*_calls.vcf.gz"), path("*_calls.vcf.gz.{tbi,csi}"), emit: calls_vcf, optional: true
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

    if (!reads_index && !ccnv) {
        error "CCNV_PREDICT: Missing reads index file"
    }

    """
    ccnv predict \\
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
    touch ${prefix}_bins.bed.gz.csi
    echo "" | gzip >  ${prefix}_calls.bed.gz
    touch ${prefix}_calls.bed.gz.tbi
    touch ${prefix}_calls.bed.gz.csi
    echo "" | gzip > ${prefix}_calls.vcf.gz
    touch ${prefix}_calls.vcf.gz.tbi
    touch ${prefix}_calls.vcf.gz.csi
    touch ${prefix}_bins.vcf.gz.tbi
    touch ${prefix}_bins.vcf.gz.csi
    touch ${prefix}_plot.html
    echo "" | gzip >  ${prefix}_segments.bed.gz
    touch ${prefix}_segments.bed.gz.tbi
    touch ${prefix}_segments.bed.gz.csi
    echo "" | gzip > ${prefix}_segments.vcf.gz
    touch ${prefix}_segments.vcf.gz.tbi
    touch ${prefix}_segments.vcf.gz.csi
    """
}
