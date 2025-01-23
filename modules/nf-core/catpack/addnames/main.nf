process CATPACK_ADDNAMES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cat:6.0.1--hdfd78af_0'
        : 'biocontainers/cat:6.0.1--hdfd78af_0'}"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(taxonomy)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${input}" == "${prefix}.txt") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    CAT_pack \\
        add_names \\
        ${args} \\
        -i ${input} \\
        -t ${taxonomy} \\
        -o ${prefix}.txt

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
        add_names \\
        ${args} \\
        -i ${input} \\
        -o ${prefix}.txt"

    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        catpack: \$(CAT_pack --version | sed 's/CAT_pack pack v//g;s/ .*//g')
    END_VERSIONS
    """
}
