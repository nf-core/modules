process FQTK {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fqtk:0.2.1--h9f5acd7_0' :
        'biocontainers/fqtk:0.2.1--h9f5acd7_0' }"

    input:
    tuple val(meta), path(sample_sheet), val(fastq_readstructure_pairs)
    // fastq_readstructure_pairs example:
    // [[<fastq name: string>, <read structure: string>, <path to fastqs: path>], [example_R1.fastq.gz, 150T, ./work/98/30bc..78y/fastqs/]]

    output:
    // Demultiplexed file name changes depending on the arg '--output-types'
    tuple val(meta), path('output/*.fq.gz')                         , emit: sample_fastq
    tuple val(meta), path('output/demux-metrics.txt')               , emit: metrics
    tuple val(meta), path('output/unmatched*.fq.gz')                , emit: most_frequent_unmatched
    path "versions.yml"                                             , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Join the absolute path from UNTAR.out.untar to the fastq file names
    fastqs = fastq_readstructure_pairs.collect{it[2]/it[0]}.join(" ")
    // Create a list of read structures, Example: 8B 8B 150T
    read_structures = fastq_readstructure_pairs.collect{it[1]}.join(" ")

    """
    mkdir output
    fqtk \\
        demux \\
            --inputs ${fastqs} \\
            --read-structures ${read_structures} \\
            --output output/ \\
            --sample-metadata ${sample_sheet} \\
            ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fqtk: \$(echo \$(fqtk --version 2>&1) | cut -d " " -f2)
    """
}
