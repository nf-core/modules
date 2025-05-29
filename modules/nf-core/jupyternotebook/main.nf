process JUPYTERNOTEBOOK {
    tag "$meta.id"
    label 'process_low'

    //NB: You likely want to override this with a container containing all required
    //dependencies for your analysis. The container at least needs to contain the
    //ipykernel, jupytext, papermill and nbconvert Python packages.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-514b1a5d280c7043110b2a8d0a87b57ba392a963:879972fc8bdc81ee92f2bce3b4805d89a772bf84-0' :
        'biocontainers/mulled-v2-514b1a5d280c7043110b2a8d0a87b57ba392a963:879972fc8bdc81ee92f2bce3b4805d89a772bf84-0' }"

    input:
    tuple val(meta), path(notebook)
    val parameters
    path input_files

    output:
    tuple val(meta), path("*.html")    , emit: report
    tuple val(meta), path("artifacts/"), emit: artifacts, optional: true
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def parametrize = (task.ext.parametrize == null) ?  true : task.ext.parametrize
    def implicit_params = (task.ext.implicit_params == null) ? true : task.ext.implicit_params
    def meta_params = (task.ext.meta_params == null) ? true : task.ext.meta_params
    def kernel   = task.ext.kernel ?: '-'

    // Dump parameters to yaml file.
    // Using a yaml file over using the CLI params because
    //  * no issue with escaping
    //  * allows to pass nested maps instead of just single values
    def params_cmd = ""
    def render_cmd = ""
    def meta_string = meta.collect { key, value -> "${key}: ${value}" }.join('\n')
    def params_string = parameters.collect { key, value -> "${key}: ${value}" }.join('\n')
    if (parametrize) {
        if (implicit_params) {
            params_cmd += 'echo cpus: ' + task.cpus + ' >> ./.params.yml \n'
            params_cmd += 'echo "artifact_dir: artifacts"  >> ./.params.yml \n'
            params_cmd += 'echo "input_dir: ./"  >> ./.params.yml \n'
        }
        if (meta_params) {
            params_cmd += 'echo "' + meta_string + '" >> ./.params.yml \n'
        }
        params_cmd += 'echo "' + params_string + '" >> ./.params.yml'
        render_cmd = "papermill -f .params.yml"
    } else {
        render_cmd = "papermill"
    }
    """
    set -o pipefail

    # write .params.yml file if required
    $params_cmd

    # Create output directory
    mkdir artifacts

    # Set parallelism for BLAS/MKL etc. to avoid over-booking of resources
    export MKL_NUM_THREADS="$task.cpus"
    export OPENBLAS_NUM_THREADS="$task.cpus"
    export OMP_NUM_THREADS="$task.cpus"
    export NUMBA_NUM_THREADS="$task.cpus"

    # Set temporary directory to remove warning about Matplotlib creating temporary directory
    export MPLCONFIGDIR=./tmp

    # Convert notebook to ipynb using jupytext, execute using papermill, convert using nbconvert
    jupytext --to notebook --output - --set-kernel ${kernel} ${notebook} > ${notebook}.ipynb
    ${render_cmd} ${notebook}.ipynb ${notebook}.executed.ipynb
    jupyter nbconvert --stdin --to html --output ${prefix}.html < ${notebook}.executed.ipynb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jupytext: \$(jupytext --version)
        ipykernel: \$(python -c "import ipykernel; print(ipykernel.__version__)")
        nbconvert: \$(jupyter nbconvert --version)
        papermill: \$(papermill --version | cut -f1 -d' ')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html
    mkdir artifacts

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jupytext: \$(jupytext --version)
        ipykernel: \$(python -c "import ipykernel; print(ipykernel.__version__)")
        nbconvert: \$(jupyter nbconvert --version)
        papermill: \$(papermill --version | cut -f1 -d' ')
    END_VERSIONS
    """
}

/**
 * Multiline code blocks need to have the same indentation level
 * as the `script:` section. This function re-indents code to the specified level.
 */
def indent_code_block(code, n_spaces) {
    def indent_str = " ".multiply(n_spaces)
    return code.stripIndent().split("\n").join("\n" + indent_str)
}

/**
 * Create a config YAML file from a groovy map
 *
 * @params task The process' `task` variable
 * @returns a line to be inserted in the bash script.
 */
def dump_params_yml(params) {
    def options = new org.yaml.snakeyaml.DumperOptions();
    options.setDefaultFlowStyle(org.yaml.snakeyaml.DumperOptions.FlowStyle.BLOCK);

    // Properly handle Groovy GStrings
    // see https://stackoverflow.com/a/35108062/2340703
    def representer = new org.yaml.snakeyaml.representer.Representer() {{
        this.multiRepresenters.put(GString, this.representers.get(String))
    }}

    def yaml = new org.yaml.snakeyaml.Yaml(representer, options)
    def yaml_str = yaml.dump(params)

    // Writing the .params.yml file directly as follows does not work.
    // It only works in 'exec:', but not if there is a `script:` section:
    // task.workDir.resolve('.params.yml').text = yaml_str

    // Therefore, we inject it into the bash script:
    return """\
        cat <<"END_PARAMS_SECTION" > ./.params.yml
        ${indent_code_block(yaml_str, 8)}
        END_PARAMS_SECTION
    """
}
