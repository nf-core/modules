process VIRALCONSENSUS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/viral_consensus:1.0.0--hcf1f8c1_0':
        'biocontainers/viral_consensus:1.0.0--hcf1f8c1_0' }"

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
    tuple val("${task.process}"), val('viralconsensus'), eval("viral_consensus --version 2>&1 | sed 's/.*v//'"), topic: versions, emit: versions_viralconsensus

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
