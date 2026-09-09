process TAXPASTA_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/taxpasta:0.7.0--pyhdfd78af_1'
        : 'quay.io/biocontainers/taxpasta:0.7.0--pyhdfd78af_1'}"

    input:
    tuple val(meta), path(profiles)
    val profiler
    val format
    path taxonomy
    path samplesheet

    output:
    tuple val(meta), path("*.{tsv,csv,arrow,parquet,biom}"), emit: merged_profiles
    tuple val("${task.process}"), val('taxpasta'), eval('taxpasta --version'), emit: versions_taxpasta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def taxonomy_option = taxonomy ? "--taxonomy ${taxonomy}" : ''
    def samplesheet_input = samplesheet ? "-s ${samplesheet}" : ''
    """
    taxpasta merge \\
        --profiler ${profiler} \\
        --output ${prefix}.${format} \\
        ${args} \\
        ${taxonomy_option} \\
        ${samplesheet_input} \\
        ${profiles}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${format}
    """
}
