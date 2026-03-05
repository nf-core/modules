process SPACERANGER_COUNT {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/spaceranger:9c5e7dc93c32448e"

    input:
    tuple val(meta), path(reads), path(image), val(slide), val(area), path(cytaimage), path(darkimage), path(colorizedimage), path(alignment), path(slidefile)
    path(reference)
    path(probeset)

    output:
    tuple val(meta), path("outs/**"), emit: outs
    tuple val("${task.process}"), val('spaceranger'), eval('spaceranger -V | sed "s/spaceranger spaceranger-//"'), emit: versions_spaceranger, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SPACERANGER_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Add flags for optional inputs on demand.
    def probeset_opt = probeset ? "--probe-set=\"${probeset}\"" : ""
    def alignment_opt = alignment ? "--loupe-alignment=\"${alignment}\"" : ""
    def slidefile_opt = slidefile ? "--slidefile=\"${slidefile}\"" : ""
    def image_opt = image ? "--image=\"${image}\"" : ""
    def cytaimage_opt = cytaimage ? "--cytaimage=\"${cytaimage}\"" : ""
    def darkimage_opt = darkimage ? "--darkimage=\"${darkimage}\"" : ""
    def colorizedimage_opt = colorizedimage ? "--colorizedimage=\"${colorizedimage}\"" : ""
    if (slide.matches("visium-(.*)") && area == "" && slidefile_opt == "") {
        slide_and_area = "--unknown-slide=\"${slide}\""
    } else {
        slide_and_area = "--slide=\"${slide}\" --area=\"${area}\""
    }
    """
    spaceranger count \\
        --id="${prefix}" \\
        --sample="${meta.id}" \\
        --fastqs=. \\
        --transcriptome="${reference}" \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        $image_opt \\
        $cytaimage_opt \\
        $darkimage_opt \\
        $colorizedimage_opt \\
        $slide_and_area \\
        $probeset_opt \\
        $alignment_opt \\
        $slidefile_opt \\
        $args
    mv ${prefix}/outs outs
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SPACERANGER_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    mkdir -p outs/
    touch outs/fake_file.txt
    """
}
