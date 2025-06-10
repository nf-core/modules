include { dumpParamsYaml ; indentCodeBlock } from "./parametrize"

// NB: You'll likely want to override this with a container containing all
// required dependencies for your analyses. Or use wave to build the container
// for you from the environment.yml You'll at least need Quarto itself,
// Papermill and whatever language you are running your analyses on; you can see
// an example in this module's environment file.
process QUARTONOTEBOOK {
    tag "${meta.id}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/28/28717ccd9ce22dbfc219f3db088d5a1fc2ca1f575b5c65621218596dcdbaac95/data'
        : 'community.wave.seqera.io/library/jupyter_matplotlib_papermill_quarto_r-rmarkdown:6d15193ce3dfc665'}"

    input:
    tuple val(meta), path(notebook)
    val parameters
    path input_files
    path extensions

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("${notebook}"), emit: notebook
    tuple val(meta), path("artifacts/*"), emit: artifacts, optional: true
    tuple val(meta), path("params.yml"), emit: params_yaml, optional: true
    tuple val(meta), path("_extensions"), emit: extensions, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def parametrize = task.ext.parametrize == null ? true : task.ext.parametrize
    def implicit_params = task.ext.implicit_params == null ? true : task.ext.implicit_params
    def meta_params = task.ext.meta_params == null ? true : task.ext.meta_params

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
    ( params_cmd ?
    """
    # Dump .params.yml heredoc (if applicable)
    ${indentCodeBlock(params_cmd, 4)}
    """ : "")
    <<
    """
    # Create output directory
    mkdir artifacts

    # Set environment variables needed for Quarto rendering
    export XDG_CACHE_HOME="./.xdg_cache_home"
    export XDG_DATA_HOME="./.xdg_data_home"

    # Fix Quarto for Apptainer (see https://community.seqera.io/t/confusion-over-why-a-tool-works-in-docker-but-fails-in-singularity-when-the-installation-doesnt-differ-i-e-using-wave-micromamba/1244)
    ENV_QUARTO=/opt/conda/etc/conda/activate.d/quarto.sh
    set +u
    if [ -z "\${QUARTO_DENO}" ] && [ -f "\${ENV_QUARTO}" ]; then
        source "\${ENV_QUARTO}"
    fi
    set -u

    # Set parallelism for BLAS/MKL etc. to avoid over-booking of resources
    export MKL_NUM_THREADS="${task.cpus}"
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export OMP_NUM_THREADS="${task.cpus}"
    export NUMBA_NUM_THREADS="${task.cpus}"

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
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Fix Quarto for Apptainer (see https://community.seqera.io/t/confusion-over-why-a-tool-works-in-docker-but-fails-in-singularity-when-the-installation-doesnt-differ-i-e-using-wave-micromamba/1244)
    ENV_QUARTO=/opt/conda/etc/conda/activate.d/quarto.sh
    set +u
    if [ -z "\${QUARTO_DENO}" ] && [ -f "\${ENV_QUARTO}" ]; then
        source "\${ENV_QUARTO}"
    fi
    set -u

    touch ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto -v)
        papermill: \$(papermill --version | cut -f1 -d' ')
    END_VERSIONS
    """
}
