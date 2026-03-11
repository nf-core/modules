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
