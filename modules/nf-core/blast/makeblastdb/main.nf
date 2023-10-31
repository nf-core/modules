process BLAST_MAKEBLASTDB {
    tag "$fasta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
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
    def is_compressed = fasta.name.endsWith(".gz")
    def fasta_name = fasta.name.replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    makeblastdb \\
        -in ${fasta_name} \\
        ${args}
    mkdir blast_db
    mv ${fasta_name}* blast_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def fasta_name = fasta.name.replace(".gz", "")
    """
    touch ${fasta_name}.fasta
    touch ${fasta_name}.fasta.ndb
    touch ${fasta_name}.fasta.nhr
    touch ${fasta_name}.fasta.nin
    touch ${fasta_name}.fasta.njs
    touch ${fasta_name}.fasta.not
    touch ${fasta_name}.fasta.nsq
    touch ${fasta_name}.fasta.ntf
    touch ${fasta_name}.fasta.nto
    mkdir blast_db
    mv ${fasta_name}* blast_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
