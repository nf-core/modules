process PLASTID_PSITE {
    tag "$meta.id"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plastid:0.6.1--py39had3e4b6_2':
        'biocontainers/plastid:0.6.1--py39had3e4b6_2' }"

    input:
    tuple val(meta), path(bam), path(bam_index)
    tuple val(meta2), path(rois_txt)

    output:
    tuple val(meta), path("*_metagene_profiles.txt"), emit: metagene_profiles
    tuple val(meta), path("*_p_offsets.png")        , emit: p_offsets_png
    tuple val(meta), path("*_p_offsets.txt")        , emit: p_offsets
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""
    """
    psite \
        "$rois_txt" \\
        "$prefix" \\
        --count_files "$bam" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plastid: 0.6.1
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_metagene_profiles.txt
    touch ${prefix}_p_offsets.png
    touch ${prefix}_p_offsets.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plastid: 0.6.1
    END_VERSIONS
    """
}
