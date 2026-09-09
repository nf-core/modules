process ESTSFS {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/est-sfs:2.04--h245ed52_0':
        'quay.io/biocontainers/est-sfs:2.04--h245ed52_0' }"

    input:
    tuple val(meta), path(e_config), path(data), path(seed)

    output:
    tuple val(meta), path("${prefix}_sfs.txt")               , emit: sfs_out
    tuple val(meta), path("${prefix}_pvalues.txt")           , emit: pvalues_out
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    tuple val("${task.process}"), val('est-sfs'), val('2.04'), emit: versions_estsfs, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    est-sfs ${e_config} ${data} ${seed} ${prefix}_sfs.txt ${prefix}_pvalues.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sfs.txt
    touch ${prefix}_pvalues.txt
    """
}
