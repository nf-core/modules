include { dumpParamsYaml; indentCodeBlock } from "./parametrize"

process QUARTONOTEBOOK {
    tag "$meta.id"
    label 'process_low'

    // NB: You'll likely want to override this with a container containing all
    // required dependencies for your analyses. You'll at least need Quarto
    // itself, Papermill and whatever language you are running your analyses on;
    // you can see an example in this module's Dockerfile.
    conda "conda-forge::jupyter=1.0.0 conda-forge::matplotlib=3.4.3 conda-forge::papermill=2.4.0 conda-forge::quarto=1.3.433 conda-forge::r-rmarkdown=2.25"
    container "docker.io/erikfas/quartonotebook"

    input:
    tuple val(meta), path(notebook)
    val parameters
    path input_files

    output:
    tuple val(meta), path("*.html")     , emit: html
    tuple val(meta), path("artifacts/*"), emit: artifacts, optional: true
    tuple val(meta), path("params.yml") , emit: params_yaml, optional: true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba on ARM64
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        arch = System.getProperty("os.arch")
        if (arch == "arm64" || arch == "aarch64") {
            exit 1, "The QUARTONOTEBOOK module does not support Conda/Mamba on ARM64. Please use Docker / Singularity / Podman instead."
        }
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def parametrize = (task.ext.parametrize == null) ?  true : task.ext.parametrize
    def implicit_params = (task.ext.implicit_params == null) ? true : task.ext.implicit_params
    def meta_params = (task.ext.meta_params == null) ? true : task.ext.meta_params

    // Dump parameters to yaml file.
    // Using a YAML file over using the CLI params because
    //  - No issue with escaping
    //  - Allows passing nested maps instead of just single values
    //  - Allows running with the language-agnostic `--execute-params`
    def params_cmd = ""
    def render_args = ""
    if (parametrize) {
        nb_params = [:]
        if (implicit_params) {
            nb_params["cpus"] = task.cpus
            nb_params["artifact_dir"] = "artifacts"
            nb_params["input_dir"] = "./"
        }
        if (meta_params) {
            nb_params["meta"] = meta
        }
        nb_params += parameters
        params_cmd = dumpParamsYaml(nb_params)
        render_args = "--execute-params params.yml"
    }
    """
    # Dump .params.yml heredoc (section will be empty if parametrization is disabled)
    ${indentCodeBlock(params_cmd, 4)}

    # Create output directory
    mkdir artifacts

    # Set parallelism for BLAS/MKL etc. to avoid over-booking of resources
    export MKL_NUM_THREADS="$task.cpus"
    export OPENBLAS_NUM_THREADS="$task.cpus"
    export OMP_NUM_THREADS="$task.cpus"
    export NUMBA_NUM_THREADS="$task.cpus"

    # Render notebook
    quarto render \\
        ${notebook} \\
        ${render_args} \\
        ${args} \\
        --output ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto -v)
        papermill: \$(papermill --version | cut -f1 -d' ')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto -v)
    END_VERSIONS
    """
}
