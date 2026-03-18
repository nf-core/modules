process COBIONTID_KMERCOUNTER {
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
    def prefix          = task.ext.prefix   ?: "${meta.id}"
    def args            = task.ext.args     ?: ""
    def is_compressed   = fasta.getExtension() == "gz" ? true : false
    def fasta_name      = is_compressed ? fasta.getBaseName() : fasta

    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    kmer-counter \\
        -f ${fasta_name} \\
        -k ${kmer_size} \\
        ${args} \\
        -o ${prefix}_kmer_counts.npy

    if [ "${is_compressed}" == "true" ]; then
        rm ${fasta_name}
    fi
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    """
    touch ${prefix}_kmer_counts.npy
    """
}
