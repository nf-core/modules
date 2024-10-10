process XENIUMRANGER_IMPORT_SEGMENTATION {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/xeniumranger:3.0.1"

    input:
    path(xenium_bundle)
    path(segmentation_file)

    output:
    path("outs/**"), emit: outs
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

    // image based segmentation
    def expansion_distance = expansion_distance ? "--expansion-distance=\"${expansion_distance}\"": ""
    def coordinate_transform = coordinate_transform ? "--coordinate-transform=\"${coordinate_transform}\"": ""
    def nuclei = nuclei ? "--nuclei=\"${segmentation_file}\"": ""
    def cells = cells ? "--cells=\"${segmentation_file}\"": ""

    // transcript based segmentation
    def transcript_assignment = transcript_assignment ? "--transcript-assignment=\"${transcript_assignment}\"": ""
    def cell_boundary_polygons = cell_boundary_polygons ? "--viz-polygons=\"${cell_boundary_polygons}\"":""

    def units = units ? "--units=\"${units}\"": ""

    """
    xeniumranger import-segmentation \\
        --id="${prefix}" \\
        --xenium-bundle="${xenium_bundle}" \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        ${nuclei} \\
        ${cells} \\
        ${expansion_distance} \\
        ${transcript_assignment} \\
        ${coordinate_transform} \\
        ${cell_boundary_polygons} \\
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
    """
    mkdir -p outs/
    touch outs/fake_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xeniumranger: \$(xeniumranger -V | sed -e "s/xeniumranger-/- /g")
    END_VERSIONS
    """
}
