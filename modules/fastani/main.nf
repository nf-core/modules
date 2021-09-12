// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FASTANI {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::fastani=1.32" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fastani:1.32--he1c1bb9_0"
    } else {
        container "quay.io/biocontainers/fastani:1.32--he1c1bb9_0"
    }

    input:
    val(meta)
    path(query)
    path(reference)
    path(query_list)
    path(reference_list)

    output:
    path("*.ani.txt")             , emit: ani
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def query_arg = query_list ? "-ql ${query_list}" : "-q ${query}"
    def reference_arg = = reference_list ? "-rl ${reference_list}" : "-r ${reference}"

    """
    fastANI \\
    ${query_arg} \\
    ${reference_arg} \\
    -o ${prefix}.ani.txt

    echo \$(fastANI --version 2>&1) | sed 's/version//;' > ${software}.version.txt
    """
}
