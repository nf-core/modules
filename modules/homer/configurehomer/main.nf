// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

def VERSION = '4.11'
def genome = 'hg19'

process HOMER_CONFIGUREHOMER {
    time '10m'
    tag "configureHomer"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:"configureHomer") }

    conda (params.enable_conda ? "bioconda::homer=4.11" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/homer:4.11--pl526h9a982cc_2"
    } else {
        container "quay.io/biocontainers/homer:4.11--pl526h9a982cc_2"
    }

    output:
    // tuple val(meta), path("/usr/local/share/homer-4.11-2/"), emit: genome
    path  "*.version.txt"              , emit: version

    script:
    def software  = getSoftwareName(task.process)
    """
    perl /usr/local/share/homer-4.11-2/configureHomer.pl \\
        -install $genome \\
        -keepScript


    echo $VERSION > ${software}.version.txt
    """
}
