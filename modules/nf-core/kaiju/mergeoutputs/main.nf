process KAIJU_MERGEOUTPUTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/kaiju:1.10.0--h43eeafb_0'
        : 'biocontainers/kaiju:1.10.0--h43eeafb_0'}"

    input:
    tuple val(meta), path(kaiju), path(kraken)
    path db

    output:
    tuple val(meta), path("*.tsv"), emit: merged
    tuple val("${task.process}"), val('kaiju'), eval("kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //'"), emit: versions_kaiju, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dbnodes = db ? '-t <(find -L ${db} -name "*nodes.dmp")' : ''

    if ("${kaiju}" == "${prefix}.tsv") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    if ("${kraken}" == "${prefix}.tsv") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    kaiju-mergeOutputs \\
        -i <(sort -k2,2 ${kaiju}) \\
        -j <(sort -k2,2 ${kraken}) \\
        -o ${prefix}.tsv \\
        ${dbnodes} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
