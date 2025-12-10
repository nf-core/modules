process PLASTID_METAGENE_GENERATE {
    tag "$annotation"
    label "process_low"

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plastid:0.6.1--py39had3e4b6_2':
        'biocontainers/plastid:0.6.1--py39had3e4b6_2' }"

    input:
    tuple val(meta), path(annotation)

    output:
    tuple val(meta), path("*_rois.txt"), emit: rois_txt
    tuple val(meta), path("*_rois.bed"), emit: rois_bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = "0.6.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    metagene generate \\
        "${annotation.baseName}" \\
        --annotation_files "$annotation" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plastid: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = "0.6.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${annotation.baseName}_rois.txt
    touch ${annotation.baseName}_rois.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plastid: $VERSION
    END_VERSIONS
    """
}
