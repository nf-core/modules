// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process DEEPTOOLS_PLOTHEATMAP {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/deeptools:3.4.3--py_0"
    //container "https://depot.galaxyproject.org/singularity/deeptools:3.4.3--py_0"

    conda (params.conda ? "bioconda::deeptools=3.4.3" : null)

    input:
    tuple val(meta), path(matrix)
    val options

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.tab"), emit: table
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    """
    plotHeatmap \\
        $ioptions.args \\
        --matrixFile $matrix \\
        --outFileName ${prefix}.plotHeatmap.pdf \\
        --outFileNameMatrix ${prefix}.plotHeatmap.mat.tab

    plotHeatmap --version | sed -e "s/plotHeatmap //g" > ${software}.version.txt
    """
}
