process OCTOPUSV_PLOTCIRCOS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/octopusv:0.4.0--pyhdfd78af_0':
        'quay.io/biocontainers/octopusv:0.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(svcf), val(output_format)
    path(fai)

    output:
    tuple val(meta), path("${prefix}.circos.${output_format}"), emit: circos
    tuple val(meta), path("${prefix}.circos.oversized_intra.tsv"), emit: oversized_intra
    tuple val("${task.process}"), val('octopusv'), eval("octopusv  --help  | sed -n 's/^ *Version: *//p'"), emit: versions_octopusv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def fai_arg = fai ? "--fai ${fai}" : ""
    def supported_formats = ['pdf', 'png', 'svg']

    if (!(output_format in supported_formats)) {
        error "Unsupported Circos output format: ${output_format}. Supported formats: ${supported_formats.join(', ')}"
    }

    """
    octopusv plot-circos \\
        --input-file ${svcf} \\
        --output-file ${prefix}.circos.${output_format} \\
        --sample ${prefix} \\
        ${fai_arg} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.circos.${output_format}
    touch ${prefix}.circos.oversized_intra.tsv
    """
}
