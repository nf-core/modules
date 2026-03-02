process MATLOCK_BAM2JUICER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/59/59df5236cc790cac47380a5654f39052a1cb3f9c7868ed397e4b3205a9fb2776/data':
        'community.wave.seqera.io/library/matlock_samtools:3c30bc2808902fde' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.links.txt")    , emit: links_txt
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '20181227'
    def args2 = task.ext.args2 ?: ''
    def filter_cmd = args2 != '' ? "samtools view --threads $task.cpus $args2 -Sb $bam > ${prefix}.juicer.bam" : ''
    def juicer_input = args2 != '' ? "${prefix}.juicer.bam" : "$bam"
    """
    $filter_cmd

    matlock \\
        bam2 \\
        juicer \\
        $juicer_input \\
        ${prefix}.links.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        matlock: $VERSION
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '20181227'
    def args2 = task.ext.args2 ?: ''
    def filter_cmd = args2 != '' ? "samtools view --threads $task.cpus $args2 -Sb $bam > ${prefix}.juicer.bam" : ''
    def juicer_input = args2 != '' ? "${prefix}.juicer.bam" : "$bam"
    """
    echo "filter command: $filter_cmd"
    echo "juicer input: $juicer_input"

    touch ${prefix}.links.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        matlock: $VERSION
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
