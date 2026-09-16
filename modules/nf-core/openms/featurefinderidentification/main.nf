process OPENMS_FEATUREFINDERIDENTIFICATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.5.0--h78fb946_0' :
        'quay.io/biocontainers/openms:3.5.0--h78fb946_0' }"

    input:
    tuple val(meta), path(mzml), path(id_int), path(id_ext), val(out_type)

    output:
    tuple val(meta), path("${prefix}.${out_type}"), emit: features
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | sed -nE 's/^Version: ([0-9.]+).*/\\1/p'"), emit: versions_openms, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def id_ext_arg = id_ext ? "-id_ext $id_ext" : ''
    prefix         = task.ext.prefix ?: "${meta.id}"
    """
    FeatureFinderIdentification \\
        -in $mzml \\
        -id $id_int \\
        $id_ext_arg \\
        -out ${prefix}.${out_type} \\
        -threads $task.cpus \\
        $args
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${out_type}
    """
}
