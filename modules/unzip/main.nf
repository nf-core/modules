// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process UNZIP {
    tag "$archive"
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
    path "${archive.baseName}/", emit: unzipped_archive
    path "versions.yml"        , emit: versions

    script:

    if ( archive instanceof List && archive.name.size > 1 ) { exit 1, "[UNZIP] error: 7za only accepts a single archive as input. Please check module input." }
    """
    7za \\
        e \\
        -o"${archive.baseName}"/ \\
        $options.args \\
        $archive

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        7za: \$(echo \$(7za --help) | sed 's/.*p7zip Version //; s/(.*//')
    END_VERSIONS
    """
}
