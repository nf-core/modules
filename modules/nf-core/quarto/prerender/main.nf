// NB 1: You'll likely want to override this with a container containing all
// required dependencies for your analyses, or use wave to build the container
// for you from the environment.yml. You'll at least need Quarto itself,
// Papermill and whatever language you are running your analyses on; you can see
// an example in this module's environment file.
//
// NB 2: You'll need to export the versions of the packages you are using inside
// your notebook to a `versions.csv` file (formatted as `package,version`),
// which will be added to the `versions` topic; module versions are handled
// separately by `eval()` statements.
process QUARTO_PRERENDER {
    tag "${prefix}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/28/28717ccd9ce22dbfc219f3db088d5a1fc2ca1f575b5c65621218596dcdbaac95/data'
        : 'community.wave.seqera.io/library/jupyter_matplotlib_papermill_quarto_r-rmarkdown:6d15193ce3dfc665'}"

    input:
    tuple val(meta), path(notebook)
    val parameters
    path input_files

    output:
    tuple val(meta), path("${prefix}{.md,_files}")                                             , emit: rendered
    tuple val(meta), path(notebook)                                                            , emit: notebook
    tuple val(meta), path("${prefix}-params.yml")                                              , emit: params_yaml
    tuple val(meta), path("${notebook_parameters.artifact_dir}/*")                             , emit: artifacts, optional: true
    path "versions.yml"                                                                        , emit: versions          , topic: versions
    tuple val("${task.process}"), val('quarto')   , eval('quarto -v')                          , emit: versions_quarto   , topic: versions
    tuple val("${task.process}"), val('papermill'), eval('papermill --version | cut -f1 -d" "'), emit: versions_papermill, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${notebook.baseName}"
    // Implicit parameters can be overwritten by supplying a value with parameters
    notebook_parameters = [
        meta: meta,
        cpus: task.cpus,
        artifact_dir: "${prefix}-artifacts",
    ] + (parameters ?: [:])
    // Parse parameters through a YAML file, which is better than CLI because:
    //  - No issue with escaping
    //  - Allows passing nested maps instead of just single values
    //  - Allows running with the language-agnostic `--execute-params`
    def yamlBuilder = new groovy.yaml.YamlBuilder()
    yamlBuilder.call(notebook_parameters)
    def yaml_content = yamlBuilder.toString().tokenize('\n').join("\n    ")
    """
    # Dump parameters to yaml file
    cat <<- END_YAML_PARAMS > ${prefix}-params.yml
    ${yaml_content}
    END_YAML_PARAMS

    # Create output directory
    mkdir "${notebook_parameters.artifact_dir}"

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

    # Render notebook to markdown
    quarto render \\
        ${notebook} \\
        ${args} \\
        --to markdown \\
        --execute-params ${prefix}-params.yml \\
        --output ${prefix}.md

    # Check that notebook package versions is exported
    if [ ! -f versions.csv ]; then
        echo "ERROR: versions.csv not found; the .qmd script must write out [tool,version] pairs used within the notebook." >&2
        exit 1
    fi

    # Write notebook package versions to YAML
    cat <<-END_VERSIONS > versions.yml
    "${task.process}_${prefix}":
    \$(awk -F',' '{printf "    %s: %s\n", \$1, \$2}' versions.csv)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${notebook.baseName}"
    notebook_parameters = [
        meta: meta,
        cpus: task.cpus,
        artifact_dir: "artifacts",
    ] + (parameters ?: [:])
    """
    # Note: The fix is also needed in the stub for `quarto -v` to work.
    ENV_QUARTO=/opt/conda/etc/conda/activate.d/quarto.sh
    set +u
    if [ -z "\${QUARTO_DENO}" ] && [ -f "\${ENV_QUARTO}" ]; then
        source "\${ENV_QUARTO}"
    fi
    set -u

    touch ${prefix}.md
    touch ${prefix}-params.yml
    touch versions.yml
    """
}
