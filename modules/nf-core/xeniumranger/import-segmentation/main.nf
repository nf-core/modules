process XENIUMRANGER_IMPORT_SEGMENTATION {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/xeniumranger:3.0.1"

    input:
    tuple val(meta), path(xenium_bundle)
    val(expansion_distance)
    path(coordinate_transform)
    path(nuclei)
    path(cells)
    path(transcript_assignment)
    path(viz_polygons)

    output:
    tuple val(meta), path("**/outs/**"), emit: outs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "XENIUMRANGER_IMPORT-SEGMENTATION module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // image based segmentation options
    def expansion_distance = expansion_distance ? "--expansion-distance=\"${expansion_distance}\"": "" // expansion distance (default - 5, range - 0 - 100)
    def coordinate_transform = coordinate_transform ? "--coordinate-transform=\"${coordinate_transform}\"": ""

    def nuclei_detection = nuclei ? "--nuclei=\"${nuclei}\"": ""
    def cells = cells ? "--cells=\"${cells}\"": ""

    // transcript based segmentation
    def transcript_assignment = transcript_assignment ? "--transcript-assignment=\"${transcript_assignment}\"": ""
    def viz_polygons = viz_polygons ? "--viz-polygons=\"${viz_polygons}\"":""

    // shared argument
    def units = coordinate_transform ? "--units=microns": "--units=pixels"

    """
    xeniumranger import-segmentation \\
        --id="${prefix}" \\
        --xenium-bundle="${xenium_bundle}" \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        ${coordinate_transform} \\
        ${nuclei_detection} \\
        ${cells} \\
        ${expansion_distance} \\
        ${transcript_assignment} \\
        ${viz_polygons} \\
        ${units} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xeniumranger: \$(xeniumranger -V | sed -e "s/xeniumranger-/- /g")
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "XENIUMRANGER_IMPORT-SEGMENTATION module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}/outs/"
    touch "${prefix}/outs/fake_file.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xeniumranger: \$(xeniumranger -V | sed -e "s/xeniumranger-/- /g")
    END_VERSIONS
    """
}
