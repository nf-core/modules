process XENIUMRANGER_IMPORT_SEGMENTATION {
    tag "${meta.id}"
    label 'process_high'

    container "nf-core/xeniumranger:4.0"

    input:
    tuple val(meta), path(xenium_bundle, stageAs: "bundle/"), path(coordinate_transform), path(nuclei), path(cells), path(transcript_assignment), path(viz_polygons), val(units)

    output:
    tuple val(meta), path("${prefix}/outs"), emit: bundle
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("XENIUMRANGER_IMPORT-SEGMENTATION module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    // image based segmentation options
    def coord_transform = coordinate_transform ? "--coordinate-transform=\"${coordinate_transform}\"" : ""
    def nuclei_detection = nuclei ? "--nuclei=\"${nuclei}\"" : ""
    def cell_detection = cells ? "--cells=\"${cells}\"" : ""

    // transcript based segmentation
    def transcript_assign = transcript_assignment ? "--transcript-assignment=\"${transcript_assignment}\"" : ""
    def polygons = viz_polygons ? "--viz-polygons=\"${viz_polygons}\"" : ""

    // shared argument
    def space = units ? "--units=${units}" : ""

    // conditional args
    def exp_dist = nuclei ? "--expansion-distance=${params.expansion_distance}" : ""

    """
    xeniumranger import-segmentation \\
        --id="${prefix}" \\
        --xenium-bundle="${xenium_bundle}" \\
        ${exp_dist} \\
        ${coord_transform} \\
        ${nuclei_detection} \\
        ${cell_detection} \\
        ${transcript_assign} \\
        ${polygons} \\
        ${space} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xeniumranger: \$(xeniumranger -V | sed -e "s/xeniumranger-/- /g")
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("XENIUMRANGER_IMPORT-SEGMENTATION module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p "${prefix}/outs"
    touch "${prefix}/outs/fake_file.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xeniumranger: \$(xeniumranger -V | sed -e "s/xeniumranger-/- /g")
    END_VERSIONS
    """
}
