process PLASTID_METAGENE_GENERATE {
    tag "$annotation"
    label "process_low"

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plastid:0.6.1--py39had3e4b6_2':
        'quay.io/biocontainers/plastid:0.6.1--py39had3e4b6_2' }"

    input:
    tuple val(meta), path(annotation)

    output:
    tuple val(meta), path("*_rois.txt"), emit: rois_txt
    tuple val(meta), path("*_rois.bed"), emit: rois_bed
    tuple val("${task.process}"), val('plastid'), val('0.6.1'), emit: versions_plastid, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    metagene generate \\
        "${annotation.baseName}" \\
        --annotation_files "$annotation" \\
        $args
    """

    stub:
    """
    touch ${annotation.baseName}_rois.txt
    touch ${annotation.baseName}_rois.bed
    """
}
