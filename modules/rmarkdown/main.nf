// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

params.options = [:]
options        = initOptions(params.options)

process RMARKDOWN {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "r-tidyverse=1.2.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/r-tidyverse:1.2.1"
    } else {
        container "quay.io/biocontainers/r-tidyverse:1.2.1"
    }

    input:
    tuple val(meta), path(notebook)
    tuple val(parameters), path(input_files)
    val(parametrize)


    output:
    tuple val(meta), path("*.html"), emit: report
    path("artifacts/*"), emit: artifacts
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    mkdir artifacts

    # export BLAS variables

    # dump parameters to yaml

    # work around  https://github.com/rstudio/rmarkdown/issues/1508
    mv "${notebook}" "${notebook}.orig"
    cp -L "${notebook}.orig" "${notebook}"

    Rscript <<EOF
    params = yaml::read_yaml(".params.yml")
    rmarkdown::render("${notebook}", params=params)
    EOF

    echo \$(Rscript -e "cat(paste(packageVersion('rmarkdown'), collapse='.'))") > ${software}.version.txt
    """
}
