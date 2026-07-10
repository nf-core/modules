process FMHFUNPROFILER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fmh-funprofiler:1.1.1--pyh106432d_0':
        'quay.io/biocontainers/fmh-funprofiler:1.1.1--pyh106432d_0' }"

    input:
    tuple val(meta), path(reads)
    path(ko_sketch)
    val(ksize)
    val(scaled)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    tuple val("${task.process}"), val('fmh-funprofiler'), eval('python -c "import importlib.metadata; print(importlib.metadata.version(\'fmh-funprofiler\'))"'), topic: versions, emit: versions_fmhfunprofiler

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    funcprofiler \\
        $args \\
        $reads \\
        $ko_sketch \\
        $ksize \\
        $scaled \\
        ${prefix}.csv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv
    """
}
