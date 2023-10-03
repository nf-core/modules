process FASTQSCREEN_FASTQSCREEN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::fastq-screen=0.15.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastq-screen%3A0.15.3--pl5321hdfd78af_0':
        'quay.io/biocontainers/fastq-screen:0.15.3--pl5321hdfd78af_0'}"

    input:
    tuple val(meta), path(reads)  // .fastq files
    path database

    output:
    tuple val(meta), path ("*_fq_screen"), emit: fastq_screen
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""

    """
    fastq_screen --threads ${task.cpus} \\
        --aligner bowtie2 \\
        --conf ${database}/fastq_screen.conf \\
        $reads \\
        $args \\
        --outdir ${prefix}_fq_screen

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqscreen: \$(echo \$(fastq_screen --version 2>&1) | sed 's/^.*FastQ Screen v//; s/ .*\$//')
    END_VERSIONS
    """
}
