process CATPACK_DOWNLOAD {
    tag "${meta.id}"
    label 'process_single'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cat:6.0.1--hdfd78af_1'
        : 'biocontainers/cat:6.0.1--hdfd78af_1'}"

    input:
    tuple val(meta), val(db)

    output:
    tuple val(meta), path("${prefix}/*.${db}.gz"), emit: fasta
    tuple val(meta), path("${prefix}/*.names.dmp"), emit: names
    tuple val(meta), path("${prefix}/*.nodes.dmp"), emit: nodes
    tuple val(meta), path("${prefix}/*accession2taxid*.gz"), emit: acc2tax
    tuple val(meta), path("${prefix}/*.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    CAT_pack \\
        download \\
        ${args} \\
        --db ${db} \\
        -o ${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        catpack: \$(CAT_pack --version | sed 's/CAT_pack pack v//g;s/ .*//g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "CAT_pack \\
        download \\
        ${args} \\
        --db ${db}
        -o ${prefix}/"

    mkdir ${prefix}/
    echo "" | gzip > ${prefix}/${prefix}.${db}.gz
    touch ${prefix}/${prefix}.names.dmp
    touch ${prefix}/${prefix}.nodes.dmp
    echo "" | gzip > ${prefix}/${prefix}.accession2taxid.gz
    touch ${prefix}/${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        catpack: \$(CAT_pack --version | sed 's/CAT_pack pack v//g;s/ .*//g')
    END_VERSIONS
    """
}
