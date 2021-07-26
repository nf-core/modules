// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
include { dump_params_yml } from "./parametrize"

params.options = [:]
options        = initOptions(params.options)
params.parametrize = true
params.implicit_params = true
params.meta_params = true

process JUPYTERNOTEBOOK {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "ipykernel jupytext nbconvert papermill" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "quay.io/biocontainers/YOUR-TOOL-HERE"
    }

    input:
    tuple val(meta), path(notebook)
    tuple val(parameters), path(input_files)


    output:
    tuple val(meta), path("*.html"), emit: report
    path("artifacts/*"), emit: artifacts
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def params_cmd = ""
    def render_cmd = ""
    if (params.parametrize) {
        nb_params = [:]
        if (params.implicit_params) {
            nb_params["cpus"] = task.cpus
            nb_params["artifact_dir"] = "artifacts"
            nb_params["input_dir"] = "."
        }
        if (params.meta_params) {
            nb_params["meta"] = meta
        }
        nb_params += parameters
        params_cmd = dump_params_yml(nb_params)
        render_cmd = "papermill -f .params.yml"
    } else {
        render_cmd = "papermill"
    }

    """
    # Create output directory
    mkdir artifacts

    # Set parallelism for BLAS/MKL etc. to avoid over-booking of resources
    export MKL_NUM_THREADS="${task.cpus}"
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export OMP_NUM_THREADS="${task.cpus}"
    export NUMBA_NUM_THREADS="${task.cpus}"

    # dump parameters to yaml
    ${params_cmd}

    # Convert notebook to ipynb using jupytext, execute using papermill, convert using nbconvert
    jupytext --to notebook --output - --set-kernel - ${notebook}  \\
        | ${render_cmd} \\
        | jupyter nbconvert --stdin --to html --output ${notebook.baseName}.html

    # TODO: show to output versions of multiple tools?
    echo \$(jupytext --version) > ${software}.version.txt
    """
}
