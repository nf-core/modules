// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

def VERSION = '4.11'

process HOMER_ANNOTATEPEAKS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::homer=4.11" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/homer:4.11--pl526h9a982cc_2"
    } else {
        container "quay.io/biocontainers/homer:4.11--pl526h9a982cc_2"
    }

    input:
    tuple val(meta), path(peak)
    path  fasta
    path  gtf

    output:
    tuple val(meta), path("*annotatePeaks.txt"), emit: txt
    path  "*.version.txt"                      , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    annotatePeaks.pl \\
        $peak \\
        $fasta \\
        $options.args \\
        -gtf $gtf \\
        -cpu $task.cpus \\
        > ${prefix}.annotatePeaks.txt

    echo $VERSION > ${software}.version.txt
    """
}
