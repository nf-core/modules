process BACKSUB {
    tag "$meta.id"
    label 'process_single'

    container "docker pull ghcr.io/schapirolabor/background_subtraction:v0.3.3"

    input:
    tuple val(meta), path(image)
    tuple val(meta2), path(markerfile)

    output:
    tuple val(meta), path("*_backsub.ome.tif"), emit: out
    tuple val(meta2), path("*_bs.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    """
    python3 /background_subtraction/background_sub.py \
        --output ${sampleName+'_backsub'}.ome.tif \
        --markerout ./markers_bs.csv \
        --raw $image \
        --marker $markerfile \
        --pixel-size
        ${Opts.moduleOpts(module, mcp)}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        backsub: \$(echo v0.3.3))
    END_VERSIONS
    """
}
