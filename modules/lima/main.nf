// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LIMA {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::lima=2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/lima:2.0.0--0"
    } else {
        container "quay.io/biocontainers/lima:2.0.0--0"
    }

    input:
    tuple val(meta), path(ccs), path(pbi)
    path primers

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "*.{bam.pbi,xml,json,clips,counts,guess,report,summary}", emit: reports
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    // def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def lima_out = ccs.toString().replaceAll(/bam$/, 'fl.bam')
    """
    lima \\
        $ccs \\
        $primers \\
        $lima_out \\
        -j $task.cpus \\
        $options.args

    lima --version | grep -e 'commit' > ${software}.version.txt
    """
}
