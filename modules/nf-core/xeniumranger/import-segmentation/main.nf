process XENIUMRANGER_IMPORT_SEGMENTATION {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/xeniumranger:4.0"

    input:
    tuple val(meta), path(xenium_bundle, stageAs: "bundle/"), path(transcript_assignment), path(viz_polygons), path(nuclei), path(cells), path(coordinate_transform), val(units)

    output:
    tuple val(meta), path("${prefix}"), emit: outs
    tuple val("${task.process}"), val("xeniumranger"), eval("xeniumranger -V | sed -e 's/.*xenium-//'"), emit: versions_xeniumranger, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "XENIUMRANGER_IMPORT_SEGMENTATION module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    // nuclei and cells are for image segmentation results
    // transcript_assignment and viz_polygons are for transcript assignment results
    // they are mutually exclusive
    if ((nuclei || cells) && (transcript_assignment || viz_polygons)) {
        error "--nuclei and --cells are for image segmentation results, which are mutually exclusive with --transcript-assignment and --viz-polygons for transcript assignment results. Please use only one of them."
    }

    def assembled_args = []
    if (task.ext.args) { assembled_args << task.ext.args.trim() }
    if (nuclei) { assembled_args << "--nuclei=\"${nuclei}\"" }
    if (cells) { assembled_args << "--cells=\"${cells}\"" }
    if (transcript_assignment) { assembled_args << "--transcript-assignment=\"${transcript_assignment}\"" }
    if (viz_polygons) { assembled_args << "--viz-polygons=\"${viz_polygons}\"" }
    if (coordinate_transform) {
        assembled_args << "--coordinate-transform=\"${coordinate_transform}\""
        // if coordinate_transform is provided, units must be microns
        assembled_args << "--units=\"microns\""
    } else if (units) {
        assembled_args << "--units=\"${units}\""
    }

    def args = assembled_args ? assembled_args.join(" \\\n        ") : ""

    """
    xeniumranger import-segmentation \\
        --id="XENIUMRANGER_IMPORT_SEGMENTATION" \\
        --xenium-bundle="${xenium_bundle}" \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        ${args}

    rm -rf "${prefix}"
    mv XENIUMRANGER_IMPORT_SEGMENTATION/outs "${prefix}"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}"
    touch "${prefix}/experiment.xenium"
    """

}
