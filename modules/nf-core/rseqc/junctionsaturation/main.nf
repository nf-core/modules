process RSEQC_JUNCTIONSATURATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:5.0.4--pyhdfd78af_1' :
        'biocontainers/rseqc:5.0.4--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  bed

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.r")  , emit: rscript
    tuple val("${task.process}"), val('rseqc'), eval('junction_saturation.py --version | sed "s/junction_saturation.py //"'), emit: versions_rseqc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    junction_saturation.py \\
        -i $bam \\
        -r $bed \\
        -o $prefix \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.junctionSaturation_plot.pdf
    touch ${prefix}.junctionSaturation_plot.r
    """
}
