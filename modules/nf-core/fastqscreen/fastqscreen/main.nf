// Adapted from https://github.com/nf-core/modules/blob/master/modules/nf-core/fastqscreen/fastqscreen/main.nf  commit 522b950
process FASTQSCREEN_FASTQSCREEN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/19/19cd241b503facdddc020294a0996a5c8b2570f19ea57b3f404618ac5988a87f/data':
        'community.wave.seqera.io/library/bowtie2_bowtie_bwa_fastq-screen_minimap2:db243c473e0c4a39'}"

    input:
    tuple val(meta), path(reads)  // .fastq files
    path database
    val aligner // 'bowtie', 'bowtie2', 'bwa', 'minimap2' or empty (defaults to minimap2)

    output:
    tuple val(meta), path("*.txt")     , emit: txt
    tuple val(meta), path("*.png")     , emit: png  , optional: true
    tuple val(meta), path("*.html")    , emit: html
    tuple val(meta), path("*.fastq.gz"), emit: fastq, optional: true
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""

    def valid_aligners = ['bowtie', 'bowtie2', 'bwa', 'minimap2']
    if (aligner && !(aligner in valid_aligners)) {
        error "Invalid aligner option: '${aligner}'. Must be one of: ${valid_aligners.join(', ')}, or empty for FastQ Screen's default (bowtie2)."
    }
    def aligner_arg = aligner ? "--aligner ${aligner}" : ''

    """
    fastq_screen --threads ${task.cpus} \\
        ${aligner_arg} \\
        --conf ${database}/fastq_screen.conf \\
        $reads \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqscreen: \$(echo \$(fastq_screen --version 2>&1) | sed 's/^.*FastQ Screen v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch test_1_screen.html
    touch test_1_screen.png
    touch test_1_screen.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqscreen: \$(echo \$(fastq_screen --version 2>&1) | sed 's/^.*FastQ Screen v//; s/ .*\$//')
    END_VERSIONS
    """

}
