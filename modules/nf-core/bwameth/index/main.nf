process BWAMETH_INDEX {
    tag "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwameth:0.2.9--pyh7e72e81_0' :
        'biocontainers/bwameth:0.2.9--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(fasta, name:"BwamethIndex/")
    val use_mem2

    output:
    tuple val(meta), path("BwamethIndex"), emit: index
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def index_cmd = use_mem2 ? "index-mem2" : "index"
    """

    bwameth.py ${index_cmd} $fasta

    rm $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwameth: \$(bwameth.py --version | cut -f2 -d" ")
    END_VERSIONS
    """

    stub:
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
