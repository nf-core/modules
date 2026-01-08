process RSEQC_SPLITBAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6f44b7933e2c2b1a340dc9485869974eb032d34e81af83716eb381964ee3e5e7/data' :
        'community.wave.seqera.io/library/rseqc_r-base:2e29d2dfda9cef15' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta2), path(bed)

    output:
    tuple val(meta), path("*.in.bam")  , emit: in_bam
    tuple val(meta), path("*.ex.bam")  , emit: ex_bam
    tuple val(meta), path("*.junk.bam"), emit: junk_bam
    tuple val("${task.process}"), val('rseqc'), eval('split_bam.py --version | sed "s/split_bam.py //"'), emit: versions_rseqc, topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    split_bam.py \\
        -i $bam \\
        -r $bed \\
        -o $prefix \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.in.bam
    touch ${prefix}.ex.bam
    touch ${prefix}.junk.bam
    """
}
