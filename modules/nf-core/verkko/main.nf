process VERKKO {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/verkko:1.4.1--h48217b1_0':
        'biocontainers/verkko:1.4.1--h48217b1_0' }"

    input:
    tuple val(meta), path(reads_pacbio)

    output:
    tuple val(meta), path("*_verkko_assembly.fasta")      , emit: assembly_fasta
    tuple val(meta), path("*_homopolymer-compressed.gfa") , emit: assembly_gfa
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    verkko \\
        $args \\
        -d verkko \\
        --hifi $reads_pacbio

    mv verkko/assembly.fasta ${prefix}_verkko_assembly.fasta
    mv verkko/assembly.homopolymer-compressed.gfa ${prefix}_homopolymer-compressed.gfa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verkko: \$(verkko --version | sed 's/verkko v//g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_verkko_assembly.fasta
    touch ${prefix}_homopolymer-compressed.gfa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verkko: \$(verkko --version | sed 's/verkko v//g')
    END_VERSIONS
    """
}
