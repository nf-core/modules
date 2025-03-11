process CATPACK_SUMMARISE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cat:6.0.1--hdfd78af_1'
        : 'biocontainers/cat:6.0.1--hdfd78af_1'}"

    input:
    tuple val(meta), path(classification)
    tuple val(meta2), path(contigs)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${classification}" == "${prefix}.txt") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    def insert_contigs = contigs ? "-c ${contigs}" : ''
    """
    CAT_pack summarise \\
        ${args} \\
        -i ${classification} \\
        ${insert_contigs} \\
        -o ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        catpack: \$(CAT_pack --version | sed 's/CAT_pack pack v//g;s/ .*//g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${classification}" == "${prefix}.txt") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    def insert_contigs = contigs ? "-c ${contigs}" : ''
    """
    echo "CAT_pack summarise \\
        ${args} \\
        -i ${classification} \\
        ${insert_contigs} \\
        -o ${prefix}.txt"

    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        catpack: \$(CAT_pack --version | sed 's/CAT_pack pack v//g;s/ .*//g')
    END_VERSIONS
    """
}
