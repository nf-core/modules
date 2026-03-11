
process PRETEXTMAP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f6b88d972aaf27c7e748e2c3b5ee86065dc5ed7824d7d13937c65844242211e2/data':
        'community.wave.seqera.io/library/htslib_pretextmap_samtools:8a29f6d0f55f98f9' }"

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
