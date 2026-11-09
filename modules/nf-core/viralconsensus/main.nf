// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VIRALCONSENSUS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/3118c77d2f8b05387921156880d31c7be0eaeaa00cc3b26d1cbc1cb8f1f2a6d5/data'
:         'community.wave.seqera.io/library/viral_consensus:1.0.0--350716b73a1e900e' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)
    path primer_bed
    val save_pos_counts
    val save_ins_counts

    output:
    tuple val(meta), path("*.consensus.fa"), emit: fasta
    tuple val(meta), path("*.pos_counts.tsv"), optional: true, emit: pos_counts
    tuple val(meta), path("*.ins_counts.json"), optional: true, emit: ins_counts
    tuple val("${task.process}"), val('viralconsensus'), eval("viral_consensus --version | sed 's/.*v//'"), topic: versions, emit: versions_viralconsensus

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def primer_arg = primer_bed ? "-p ${primer_bed}" : ''
    def pos_counts_arg = save_pos_counts ? "-op ${prefix}.pos_counts.tsv" : ''
    def ins_counts_arg = save_ins_counts ? "-oi ${prefix}.ins_counts.json" : ''
    """
    viral_consensus \\
        -i $bam \\
        -r $fasta \\
        -o ${prefix}.consensus.fa \\
        $primer_arg \\
        $pos_counts_arg \\
        $ins_counts_arg \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def touch_pos_counts = save_pos_counts ? "touch ${prefix}.pos_counts.tsv" : ''
    def touch_ins_counts = save_ins_counts ? "touch ${prefix}.ins_counts.json" : ''
    """
    touch ${prefix}.consensus.fa
    $touch_pos_counts
    $touch_ins_counts
    """
}
