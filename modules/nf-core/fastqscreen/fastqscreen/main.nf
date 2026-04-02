process FASTQSCREEN_FASTQSCREEN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fc/fc53eee7ca23c32220a9662fbb63c67769756544b6d74a1ee85cf439ea79a7ee/data'
        : 'community.wave.seqera.io/library/fastq-screen_perl-gdgraph:5c1786a5d5bc1309'}"

    input:
    tuple val(meta), path(reads)
    path database

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*.png"), emit: png, optional: true
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.fastq.gz"), emit: fastq, optional: true
    tuple val("${task.process}"), val('fastqscreen'), eval('fastq_screen --version 2>&1 | sed "s/^.*FastQ Screen v//;"'), emit: versions_fastqscreen, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""

    """
    fastq_screen --threads ${task.cpus} \\
        --conf ${database}/fastq_screen.conf \\
        ${reads} \\
        ${args} \\
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}_screen.html
    touch ${prefix}_screen.png
    touch ${prefix}_screen.txt
    """
}
