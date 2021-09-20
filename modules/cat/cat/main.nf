// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CAT_CAT {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::pigz=2.3.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pigz:2.3.4"
    } else {
        container "quay.io/biocontainers/pigz:2.3.4"
    }

    input:
    path files

    output:
    path "file*"        , emit: file
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    cpus = Math.floor(task.cpus/2).toInteger()

    // Use options.suffix if specified, otherwise .out; add .gz if first input file has it
    suffix  = options.suffix ? "${options.suffix}" : ".out"
    suffix += files[0].name =~ /\.gz/ ? '.gz' : ''

    """
    cat ${options.args} $files ${options.args2} > file${suffix}
    cat --version | grep 'GNU coreutils' | sed 's/cat (GNU coreutils) //' > ${software}.version.txt
    """
}
