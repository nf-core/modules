process BLAST_MAKEBLASTDB {
    tag "$fasta"
    label 'process_medium'

    conda "bioconda::blast=2.14.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.14.1--pl5321h6f7f691_0':
        'biocontainers/blast:2.14.1--pl5321h6f7f691_0' }"

    input:
    path fasta

    output:
    path 'blast_db'    , emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    makeblastdb \\
        -in $fasta \\
        $args
    mkdir blast_db
    mv ${fasta}* blast_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch ${fasta}.fasta
    touch ${fasta}.fasta.ndb
    touch ${fasta}.fasta.nhr
    touch ${fasta}.fasta.nin
    touch ${fasta}.fasta.njs
    touch ${fasta}.fasta.not
    touch ${fasta}.fasta.nsq
    touch ${fasta}.fasta.ntf
    touch ${fasta}.fasta.nto
    mkdir blast_db
    mv ${fasta}* blast_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
