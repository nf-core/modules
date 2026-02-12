process BLAST_MAKEBLASTDB {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0c/0c86cbb145786bf5c24ea7fb13448da5f7d5cd124fd4403c1da5bc8fc60c2588/data':
        'community.wave.seqera.io/library/blast:2.17.0--d4fb881691596759' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}"), emit: db
    tuple val("${task.process}"), val("makeblastdb"), eval("makeblastdb -version 2>&1 | sed 's/^.*makeblastdb: //; s/ .*\$//'"), topic: versions, emit: versions_makeblastdb

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

    """

    stub:
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

    """
}
