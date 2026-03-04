process SHAPEIT5_PHASECOMMON {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/shapeit5:5.1.1--hb60d31d_0'
        : 'biocontainers/shapeit5:5.1.1--hb60d31d_0'}"

    input:
    tuple val(meta), path(input), path(input_index), path(pedigree), val(region), path(reference), path(reference_index), path(scaffold), path(scaffold_index), path(map)

    output:
    tuple val(meta), path("*.{bcf,graph,bh}"), emit: phased_variant
    tuple val("${task.process}"), val('shapeit5'), eval('SHAPEIT5_phase_common | sed "5!d;s/^.*Version *: //; s/ .*$//"'), topic: versions, emit: versions_shapeit5

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def extension = args.contains("--output-format bcf")   ? "bcf"   :
                    args.contains("--output-format graph") ? "graph" :
                    args.contains("--output-format bh")    ? "bh"    :
                    "bcf"

    if ("${input}" == "${prefix}.${extension}") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }

    def map_command       = map       ? "--map ${map}"             : ""
    def reference_command = reference ? "--reference ${reference}" : ""
    def scaffold_command  = scaffold  ? "--scaffold ${scaffold}"   : ""
    def pedigree_command  = pedigree  ? "--pedigree ${pedigree}"   : ""

    """
    SHAPEIT5_phase_common \\
        ${args} \\
        --input ${input} \\
        ${map_command} \\
        ${reference_command} \\
        ${scaffold_command} \\
        ${pedigree_command} \\
        --region ${region} \\
        --thread ${task.cpus} \\
        --output ${prefix}.${extension}
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def extension = args.contains("--output-format bcf")   ? "bcf"   :
                    args.contains("--output-format graph") ? "graph" :
                    args.contains("--output-format bh")    ? "bh"    :
                    "bcf"
    """
    touch ${prefix}.${extension}
    """
}
