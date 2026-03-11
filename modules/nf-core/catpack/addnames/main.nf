process CATPACK_ADDNAMES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cat:6.0.1--hdfd78af_1'
        : 'biocontainers/cat:6.0.1--hdfd78af_1'}"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(taxonomy)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: txt
    tuple val("${task.process}"), val('catpack'), eval("CAT_pack --version | sed 's/CAT_pack pack v//g;s/ .*//g'"), topic: versions, emit: versions_catpack

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${input}" == "${prefix}.txt") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    CAT_pack add_names \\
        ${args} \\
        -i ${input} \\
        -t ${taxonomy} \\
        -o ${prefix}.txt
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "CAT_pack add_names \\
        ${args} \\
        -i ${input} \\
        -o ${prefix}.txt"

    touch ${prefix}.txt
    """
}
