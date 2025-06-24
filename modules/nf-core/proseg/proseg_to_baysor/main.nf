process PROSEG_TO_BAYSOR {
    tag "$meta.id"
    label 'process_low'

    // if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
    //     error "proseg module does not support Conda. Please use Docker / Singularity / Podman instead."
    // }
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/ruijintracyyang/proseg:v1.0':
        'docker.io/ruijintracyyang/proseg:v1.0' }"

    input:
    tuple val(meta), path(transcript_metadata)
    tuple val(meta), path(cell_polygons)

    output:
    tuple val(meta), path("*baysor-cell-polygons.geojson")  , emit: baysor_cell_polygons
    tuple val(meta), path("*baysor-transcript-metadata.csv"), emit: baysor_transcript_metadata
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    proseg-to-baysor \\
        ${transcript_metadata} \\
        ${cell_polygons} \\
        --output-transcript-metadata ${prefix}-baysor-transcript-metadata.csv \\
        --output-cell-polygons ${prefix}-baysor-cell-polygons.geojson \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proseg: \$(proseg --version | sed 's/proseg //')
    END_VERSIONS
    """

    stub:
    """
    touch proseg-to-baysor-transcript-metadata.csv
    touch proseg-to-baysor-cell-polygons.geojson

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proseg: \$(proseg --version | sed 's/proseg //')
    END_VERSIONS
    """
}
