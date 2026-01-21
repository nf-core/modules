
process PRETEXTMAP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f3591ce8609c7b3b33e5715333200aa5c163aa61%3A44321ab4d64f0b6d0c93abbd1406369d1b3da684':
        'community.wave.seqera.io/library/pretextmap_samtools_pip_htlib:ea913e7946e65c1a' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.pretext")  , emit: pretext
    tuple val("${task.process}"), val('PretextMap'), eval('PretextMap | sed "/Version/!d; s/.*Version //"'), emit: versions_pretextmap, topic: versions
    tuple val("${task.process}"), val('samtools'), eval('samtools --version | sed "1!d; s/samtools //"'), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args     ?: '' // PretextMap args
    def args2         = task.ext.args2    ?: '' // Samtools view args
    def prefix        = task.ext.prefix   ?: "${meta.id}"
    def reference     = fasta             ? "--reference ${fasta}" : ""
    def pairs_input   = input.toString().endsWith(".pairs.gz")
    def input_command = pairs_input ? "zcat ${input}" : "samtools view $args2 $reference -h ${input}"
    """
    ${input_command} | PretextMap \\
        $args \\
        -o ${prefix}.pretext
    """

    stub:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pretext
    """
}
