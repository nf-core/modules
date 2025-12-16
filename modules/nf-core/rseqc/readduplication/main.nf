process RSEQC_READDUPLICATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:5.0.4--pyhdfd78af_1' :
        'biocontainers/rseqc:5.0.4--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*seq.DupRate.xls"), emit: seq_xls
    tuple val(meta), path("*pos.DupRate.xls"), emit: pos_xls
    tuple val(meta), path("*.pdf")           , emit: pdf
    tuple val(meta), path("*.r")             , emit: rscript
    tuple val("${task.process}"), val('rseqc'), eval('read_duplication.py --version | sed "s/read_duplication.py //"'), emit: versions_rseqc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    read_duplication.py \\
        -i $bam \\
        -o $prefix \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.seq.DupRate.xls
    touch ${prefix}.pos.DupRate.xls
    touch ${prefix}.DupRate_plot.pdf
    touch ${prefix}.DupRate_plot.r
    """
}
