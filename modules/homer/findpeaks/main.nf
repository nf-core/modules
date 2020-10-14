// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

def VERSION = '4.11'

process HOMER_MAKETAGDIRECTORY {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/homer:4.11--pl526h9a982cc_2"
    //container "https://depot.galaxyproject.org/singularity/homer:4.11--pl526h9a982cc_2"

    conda (params.conda ? "bioconda::homer=4.11" : null)

    input:
    tuple val(meta), path(tagDir)
    val   options

    output:
    tuple val(meta), path("*peaks.txt"), emit: txt
    path  "*.version.txt"            , emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    """
    findPeaks \\
        $tagDir \\
        -style $ioptions.style \\
        $ioptions.args \\
        -cpu $task.cpus \\
        > ${prefix}.peaks.txt

    echo $VERSION > ${software}.version.txt
    """
}
