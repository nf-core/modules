process RSEQC_READDISTRIBUTION {
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
    tuple val(meta), path("*.read_distribution.txt"), emit: txt
    tuple val("${task.process}"), val('rseqc'), eval('read_distribution.py --version | sed "s/read_distribution.py //"'), emit: versions_rseqc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    read_distribution.py \\
        -i $bam \\
        -r $bed \\
        > ${prefix}.read_distribution.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.read_distribution.txt
    """
}
