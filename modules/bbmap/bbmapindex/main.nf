// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '38.92'       // There seems to be no way to get the version out of the bbmap.sh script, only a line with "Last modified" and a date

process BBMAP_BBMAPINDEX {
    tag "$fasta"
    label 'process_long'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bbmap=${VERSION}" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bbmap:${VERSION}--he522d1c_0"
    } else {
        container "quay.io/biocontainers/bbmap:${VERSION}--he522d1c_0"
    }

    input:
    tuple val(meta), path(fasta)

    output:
    path 'ref'                    , emit: db
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    bbmap.sh \\
        ref=${fasta} \\
        $options.args \\
        threads=$task.cpus \\
        -Xmx${task.memory.toGiga()}g

    echo ${VERSION} > ${software}.version.txt
    """
}
