process SDRFPIPELINES_CONVERTMHCQUANT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sdrf-pipelines:0.1.3--pyhdfd78af_0':
        'quay.io/biocontainers/sdrf-pipelines:0.1.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(sdrf)

    output:
    tuple val(meta), path("${prefix}_samplesheet.tsv")   , emit: samplesheet
    tuple val(meta), path("${prefix}_search_presets.tsv"), emit: search_presets
    tuple val("${task.process}"), val('sdrf-pipelines'), eval("parse_sdrf --version | cut -d ' ' -f 2"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_sdrf convert-mhcquant \\
        $args \\
        -s ${sdrf} \\
        -os ${prefix}_samplesheet.tsv \\
        -op ${prefix}_search_presets.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_samplesheet.tsv
    touch ${prefix}_search_presets.tsv
    """
}
