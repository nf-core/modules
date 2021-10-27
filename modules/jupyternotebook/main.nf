// Import generic module functions
include { initOptions; saveFiles; getProcessName; getSoftwareName } from './functions'
include { dump_params_yml; indent_code_block } from "./parametrize"

params.options         = [:]
options                = initOptions(params.options)
params.parametrize     = true
params.implicit_params = true
params.meta_params     = true

process JUPYTERNOTEBOOK {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    //NB: You likely want to override this with a container containing all required
    //dependencies for your analysis. The container at least needs to contain the
    //ipykernel, jupytext, papermill and nbconvert Python packages.
    conda (params.enable_conda ? "ipykernel=6.0.3 jupytext=1.11.4 nbconvert=6.1.0 papermill=2.3.3 matplotlib=3.4.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-514b1a5d280c7043110b2a8d0a87b57ba392a963%3A879972fc8bdc81ee92f2bce3b4805d89a772bf84-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-514b1a5d280c7043110b2a8d0a87b57ba392a963:879972fc8bdc81ee92f2bce3b4805d89a772bf84-0"
    }

    input:
    tuple val(meta), path(notebook)
    val parameters
    path input_files

    output:
    tuple val(meta), path("*.html"), emit: report
    tuple val(meta), path("artifacts/"), emit: artifacts, optional: true
    path "versions.yml"            , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    // Dump parameters to yaml file.
    // Using a yaml file over using the CLI params because
    //  * no issue with escaping
    //  * allows to pass nested maps instead of just single values
    def params_cmd = ""
    def render_cmd = ""
    if (params.parametrize) {
        nb_params = [:]
        if (params.implicit_params) {
            nb_params["cpus"] = task.cpus
            nb_params["artifact_dir"] = "artifacts"
            nb_params["input_dir"] = "./"
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
    set -o pipefail

    # Dump .params.yml heredoc (section will be empty if parametrization is disabled)
    ${indent_code_block(params_cmd, 4)}

    # Create output directory
    mkdir artifacts

    # Set parallelism for BLAS/MKL etc. to avoid over-booking of resources
    export MKL_NUM_THREADS="${task.cpus}"
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export OMP_NUM_THREADS="${task.cpus}"
    export NUMBA_NUM_THREADS="${task.cpus}"

    # Convert notebook to ipynb using jupytext, execute using papermill, convert using nbconvert
    jupytext --to notebook --output - --set-kernel - ${notebook}  \\
        | ${render_cmd} \\
        | jupyter nbconvert --stdin --to html --output ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        jupytext: \$(jupytext --version)
        ipykernel: \$(python -c "import ipykernel; print(ipykernel.__version__)")
        nbconvert: \$(jupyter nbconvert --version)
        papermill: \$(papermill --version | cut -f1 -d' ')
    END_VERSIONS
    """
}
