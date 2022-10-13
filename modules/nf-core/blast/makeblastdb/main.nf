process BLAST_MAKEBLASTDB {
    tag "$fasta"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::blast=2.12.0' : null)
    def container_image = "/blast:2.12.0--pl5262h3289130_0"
                                                    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    path fasta

    output:
    path 'blast_db'     , emit: db
    path "versions.yml" , emit: versions

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
}
