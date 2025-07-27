process BLAST_MAKEBLASTDB {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_5':
        'biocontainers/blast:2.16.0--h66d330f_5' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}"), emit: db
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args ?: ''
    prefix             = task.ext.prefix ?: "${meta.id}"
    def is_compressed  = fasta.getExtension() == "gz" ? true : false
    def fasta_name     = is_compressed ? fasta.getBaseName() : fasta
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    makeblastdb \\
        -in ${fasta_name} \\
        -out ${prefix}/${fasta_name} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(makeblastdb -version 2>&1 | sed 's/^.*makeblastdb: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args           = task.ext.args ?: ''
    prefix             = task.ext.prefix ?: "${meta.id}"
    def is_compressed  = fasta.getExtension() == "gz" ? true : false
    def fasta_name     = is_compressed ? fasta.getBaseName() : fasta
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
    mkdir ${prefix}
    mv ${fasta_name}* ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(makeblastdb -version 2>&1 | sed 's/^.*makeblastdb: //; s/ .*\$//')
    END_VERSIONS
    """
}
