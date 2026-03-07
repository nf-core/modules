
process PRETEXTMAP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e9e6a49d9810ef0101a4a003afeda9b32c1d0d06b196ec13a5c9f5919bd1869e/data':
        'community.wave.seqera.io/library/htslib_pretextmap_samtools:6d973e19ac7b0a1f' }"

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
