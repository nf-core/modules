process SDRFPIPELINES_PARSESDRF {
    tag "${sdrf.baseName}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sdrf-pipelines:0.1.3--pyhdfd78af_0':
        'quay.io/biocontainers/sdrf-pipelines:0.1.3--pyhdfd78af_0' }"

    input:
    path sdrf

    output:
    path "samplesheet.tsv"   , emit: samplesheet
    path "search_presets.tsv", emit: search_presets
    tuple val("${task.process}"), val('sdrf-pipelines'), eval("parse_sdrf --version | cut -d ' ' -f 2"), topic: versions, emit: versions_sdrfpipelines

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    parse_sdrf $args \\
        -s ${sdrf}
    """

    stub:
    """
    touch samplesheet.tsv
    touch search_presets.tsv
    """
}
