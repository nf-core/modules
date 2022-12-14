process GENMAP_INDEX {
    tag "$fasta"
    label 'process_high'

    conda "bioconda::genmap=1.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genmap:1.3.0--h1b792b2_1' :
        'quay.io/biocontainers/genmap:1.3.0--h1b792b2_1' }"

    input:
    path fasta

    output:
    path "genmap"       , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    genmap \\
        index \\
        -F $fasta \\
        -I genmap

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmap: \$(genmap --version 2>&1 | sed 's/GenMap version: //; s/SeqAn.*\$//')
    END_VERSIONS
    """
}
