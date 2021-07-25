// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
include { dump_params_yml } from "./parametrize"

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

    //NB: You likely want to override this with a container containing all required
    //dependencies for you analysis. The container at least needs to contain the
    //yaml and rmarkdown R packages.
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
    val(implicit_params)
    val(meta_params)


    output:
    tuple val(meta), path("*.html"), emit: report
    path("artifacts/*"), emit: artifacts, optional: true
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def parametrize = params.parametrize ?: true
    def implicit_params = params.implicit_params ?: true
    def meta_params = params.meta_params ?: true

    def params_cmd = ""
    def render_cmd = ""
    if (parametrize) {
        nb_params = [:]
        if (implicit_params) {
            nb_params["cpus"] = task.cpus
            nb_params["artifact_dir"] = "artifacts"
            nb_params["input_dir"] = "."
        }
        if (meta_params) {
            nb_params["meta"] = meta
        }
        nb_params += parameters
        params_cmd = dump_params_yml(nb_params)
        render_cmd = (
            "params = yaml::read_yaml('.params.yml')\n" +
            "rmarkdown::render('${notebook}', params=params)"
        )
    } else {
        render_cmd = "rmarkdown::render('${notebook}')"
    }


    """
    # Create output directory
    mkdir artifacts

    # Set parallelism for BLAS/MKL etc. to avoid over-booking of resources
    export MKL_NUM_THREADS="${task.cpus}"
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export OMP_NUM_THREADS="${task.cpus}"

    # dump parameters to yaml
    ${params_cmd}

    # work around  https://github.com/rstudio/rmarkdown/issues/1508
    mv "${notebook}" "${notebook}.orig"
    cp -L "${notebook}.orig" "${notebook}"

    # Render notebook

    Rscript <<EOF
    ${render_cmd}
    EOF

    echo \$(Rscript -e "cat(paste(packageVersion('rmarkdown'), collapse='.'))") > ${software}.version.txt
    """
}
