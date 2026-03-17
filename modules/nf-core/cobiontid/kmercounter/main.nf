process KMERCOUNTER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kmer-counter:0.1.2--h4349ce8_0' :
        'biocontainers/kmer-counter:0.1.2--h4349ce8_0' }"

    input:
    tuple val(meta), path(fasta)
    val kmer_size

    output:
    tuple val(meta), path( "*.npy" ) , emit: npy
    tuple val("${task.process}"), val('kmer-counter'), eval('kmer-counter --version | sed -e "s/K-mer counter //g"'), emit: versions_kmercounter, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    def args    = task.ext.args     ?: ""
    """
    kmer-counter \\
        -f ${fasta} \\
        -k ${kmer_size} \\
        ${args} \\
        -o ${prefix}_kmer_counts.npy
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    """
    touch ${prefix}_kmer_counts.npy
    """
}
