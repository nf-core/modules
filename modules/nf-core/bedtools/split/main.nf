process BEDTOOLS_SPLIT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0'
        : 'biocontainers/bedtools:2.31.1--hf5e1c6e_0'}"

    input:
    tuple val(meta), path(bed), val(count)

    output:
    tuple val(meta), path("*.bed"), emit: beds
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), topic: versions, emit: versions_bedtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedtools \\
        split \\
        ${args} \\
        -n ${count} \\
        -i ${bed} \\
        -p ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def create_beds = (1..count)
        .collect { number ->
            def numberString = "0".multiply(4 - number.toString().size()) + "${number}"
            "    touch ${prefix}.${numberString}.bed"
        }
        .join("\n")
    """
    ${create_beds}
    """
}
