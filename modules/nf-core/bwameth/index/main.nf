process BWAMETH_INDEX {
    tag "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwameth:0.2.7--pyh7cba7a3_0' :
        'biocontainers/bwameth:0.2.7--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(fasta, name:"BwamethIndex/")

    output:
    tuple val(meta), path("BwamethIndex"), emit: index
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """

    bwameth.py index $fasta

    rm $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwameth: \$(bwameth.py --version | cut -f2 -d" ")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    rm $fasta

    mkdir -p BwamethIndex/
    touch BwamethIndex/genome.fasta.bwameth.c2t
    touch BwamethIndex/genome.fasta.bwameth.c2t.amb
    touch BwamethIndex/genome.fasta.bwameth.c2t.ann
    touch BwamethIndex/genome.fasta.bwameth.c2t.bwt
    touch BwamethIndex/genome.fasta.bwameth.c2t.pac
    touch BwamethIndex/genome.fasta.bwameth.c2t.sa


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwameth: \$(bwameth.py --version | cut -f2 -d" ")
    END_VERSIONS
    """
}
