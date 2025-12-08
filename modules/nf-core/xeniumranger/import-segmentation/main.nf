process XENIUMRANGER_IMPORT_SEGMENTATION {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/xeniumranger:4.0"

    input:
    tuple val(meta), path(xenium_bundle, stageAs: "bundle/"), path(transcript_assignment), path(viz_polygons), path(nuclei), path(cells), path(coordinate_transform)

    output:
    tuple val(meta), path("${prefix}"), emit: outs
    tuple val("${task.process}"), val("xeniumranger"), eval("xeniumranger -V | sed -e 's/xeniumranger-/- /g'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "XENIUMRANGER_IMPORT_SEGMENTATION module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    prefix = task.ext.prefix ?: "${meta.id}"
    def assembled_args = []
    if (task.ext.args) { assembled_args << task.ext.args.trim() }
    if (nuclei) { assembled_args << "--nuclei=\"${nuclei}\"" }
    if (cells) { assembled_args << "--cells=\"${cells}\"" }
    if (transcript_assignment) { assembled_args << "--transcript-assignment=\"${transcript_assignment}\"" }
    if (viz_polygons) { assembled_args << "--viz-polygons=\"${viz_polygons}\"" }
    if (coordinate_transform) {
        assembled_args << "--coordinate-transform=\"${coordinate_transform}\""
        assembled_args = assembled_args.replaceAll("--units=pixels", "--units=microns")
        }
    def args = assembled_args ? assembled_args.join(" \\\n        ") : ""
    
    """
    xeniumranger import-segmentation \\
        --id="${prefix}" \\
        --xenium-bundle="${xenium_bundle}" \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        ${args}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "XENIUMRANGER_IMPORT_SEGMENTATION module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    def assembled_args = []
    if (task.ext.args) { assembled_args << task.ext.args.trim() }
    if (nuclei) { assembled_args << "--nuclei=\"${nuclei}\"" }
    if (cells) { assembled_args << "--cells=\"${cells}\"" }
    if (transcript_assignment) { assembled_args << "--transcript-assignment=\"${transcript_assignment}\"" }
    if (viz_polygons) { assembled_args << "--viz-polygons=\"${viz_polygons}\"" }
    if (coordinate_transform) { assembled_args << "--coordinate-transform=\"${coordinate_transform}\"" }
    def args = assembled_args.join(" \\\n        ")
    
    """
    xeniumranger import-segmentation \\
        --id="${prefix}" \\
        --xenium-bundle="${xenium_bundle}" \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        ${args} \\
        --dry

    if [ -d "XENIUMRANGER_IMPORT_SEGMENTATION/outs" ]; then
        rm -rf "${prefix}"
        mv XENIUMRANGER_IMPORT_SEGMENTATION/outs "${prefix}"
    else
        mkdir -p "${prefix}"
        touch "${prefix}/dry_run.txt"
    fi
    """
}
