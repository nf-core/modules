// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process UNZIP {
    tag '$bam'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }


    conda (params.enable_conda ? "bioconda::p7zip=15.09" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/p7zip:15.09--h2d50403_4"
    } else {
        container "quay.io/biocontainers/p7zip:15.09--h2d50403_4"
    }

    input:
    path archive

    output:
    path "${archive.baseName}/" , emit: unzipped_archive
    path "*.version.txt"     , emit: version

    script:
    def software = getSoftwareName(task.process)

    if ( archive instanceof List && archive.name.size > 1 ) { exit 1, "[UNZIP] error: 7za only accepts a single archive as input. Please check module input." }

    """
    7za \\
        e \\
        -o"${archive.baseName}"/ \\
        $options.args \\
        $archive

    echo \$(7za --help) | grep Version | sed 's/.*p7zip Version//; s/(.*//' 1> ${software}.version.txt
    """
}
