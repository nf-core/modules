// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BBMAP_INDEX {
    tag "$fasta"
    label 'process_long'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::bbmap=38.92" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bbmap:38.92--he522d1c_0"
    } else {
        container "quay.io/biocontainers/bbmap:38.92--he522d1c_0"
    }

    input:
    path fasta

    output:
    path 'ref'                    , emit: index
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    bbmap.sh \\
        ref=${fasta} \\
        $options.args \\
        threads=$task.cpus \\
        -Xmx${task.memory.toGiga()}g

    echo \$(bbversion.sh) > ${software}.version.txt
    """
}
