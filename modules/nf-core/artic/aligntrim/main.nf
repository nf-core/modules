// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ARTIC_ALIGNTRIM {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/23/23fb96f7b8c846ee7e502520801ab924560be07d9c5eec04e0c804d0c44851e1/data'
:         'community.wave.seqera.io/library/samtools_pip_align_trim:5fd2b1772c911075' }"

    input:
    tuple val(meta), path(samfile), path(scheme_bed), val(normalise_depth)
    val sort_bam

    output:
    tuple val(meta), path("*.primertrimmed*.bam"),              emit: primertrimmed_bam
    tuple val(meta), path("*.align_trim_report.tsv"),           emit: align_trim_report
    tuple val(meta), path("*.amp_depth_report.tsv"),            emit: amp_depth_report
    tuple val(meta), path("*.pre-normalisation.coverage.tsv"),  emit: pre_normalisation_coverage_report
    tuple val(meta), path("*.post-normalisation.coverage.tsv"), emit: post_normalisation_coverage_report, optional: true
    tuple val(meta), path("*.summary.tsv"),                     emit: qc_summary
    tuple val(meta), path("*.amplicons.tsv"),                   emit: qc_amplicons
    tuple val(meta), path("*.primer_pairs.tsv"),                emit: qc_primer_pairs
    tuple val(meta), path("*.split_reads.tsv"),                 emit: qc_split_reads
    tuple val(meta), path("*.overlaps.tsv"),                    emit: qc_overlaps
    tuple val(meta), path("*.rejected.bam"),                    emit: rejected_bam
    tuple val("${task.process}"), val('align_trim'), eval("align_trim --version | sed 's/align_trim //'"), emit: versions_align_trim, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def normalise = normalise_depth ? "--normalise ${normalise_depth}" : ''
    def sort_bam_cmd = sort_bam ? "samtools sort -o ${prefix}.primertrimmed.sorted.bam ${prefix}.primertrimmed.bam && rm ${prefix}.primertrimmed.bam" : ''
    """
    align_trim \\
        --samfile ${samfile} \\
        --output ${prefix}.primertrimmed.bam \\
        --report ${prefix}.align_trim_report.tsv \\
        --amp-depth-report ${prefix}.amp_depth_report.tsv \\
        --qc-prefix ${prefix} \\
        ${normalise} \\
        --rejected-output ${prefix}.rejected.bam \\
        ${args} \\
        ${scheme_bed}

    ${sort_bam_cmd}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sort_bam_cmd = sort_bam ? "touch ${prefix}.primertrimmed.sorted.bam" : "touch ${prefix}.primertrimmed.bam"
    def post_normalisation_cmd = normalise_depth ? "touch ${prefix}.post-normalisation.coverage.tsv" : ''
    """
    echo ${args}
    ${sort_bam_cmd}
    touch ${prefix}.rejected.bam
    ${post_normalisation_cmd}
    touch ${prefix}.align_trim_report.tsv
    touch ${prefix}.amp_depth_report.tsv
    touch ${prefix}.pre-normalisation.coverage.tsv
    touch ${prefix}.summary.tsv
    touch ${prefix}.amplicons.tsv
    touch ${prefix}.primer_pairs.tsv
    touch ${prefix}.split_reads.tsv
    touch ${prefix}.overlaps.tsv
    """
}
