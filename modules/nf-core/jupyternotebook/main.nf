process JUPYTERNOTEBOOK {
    tag "$meta.id"
    label 'process_low'

    //NB: You likely want to override this with a container containing all required
    //dependencies for your analysis. The container at least needs to contain the
    //ipykernel, jupytext, papermill and nbconvert Python packages.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-514b1a5d280c7043110b2a8d0a87b57ba392a963:879972fc8bdc81ee92f2bce3b4805d89a772bf84-0' :
        'quay.io/biocontainers/mulled-v2-514b1a5d280c7043110b2a8d0a87b57ba392a963:879972fc8bdc81ee92f2bce3b4805d89a772bf84-0' }"

    input:
    tuple val(meta), path(notebook)
    val parameters
    path input_files
    val kernel_

    output:
    tuple val(meta), path("*.html")    , emit: report
    tuple val(meta), path("artifacts/"), emit: artifacts, optional: true
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def kernel   = kernel_ ?: '-'
    // Implicit parameters can be overwritten by supplying a value with parameters
    notebook_parameters = [
        meta: meta,
        cpus: task.cpus,
    ] + (parameters ?: [:])

    // Dump parameters to yaml file.
    // Using a yaml file over using the CLI params because
    //  * no issue with escaping
    //  * allows to pass nested maps instead of just single values
    def yamlBuilder = new groovy.yaml.YamlBuilder()
    yamlBuilder.call(notebook_parameters)
    def yaml_content = yamlBuilder.toString().tokenize('\n').join("\n    ")
    """
    # Dump parameters to yaml file
    cat <<- END_YAML_PARAMS > params.yml
    ${yaml_content}
    END_YAML_PARAMS

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
    papermill -f params.yml ${notebook}.ipynb ${notebook}.executed.ipynb
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
