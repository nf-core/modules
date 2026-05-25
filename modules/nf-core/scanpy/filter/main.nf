process SCANPY_FILTER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/8158c8824afb7a57da5327fbd188082d24d205e078103ca249d74e93cc1cd603/data'
        : 'community.wave.seqera.io/library/python_pyyaml_scanpy:da8fc259e2b95ada'}"

    input:
    tuple val(meta), path(anndata)
    val min_genes
    val min_cells
    val min_counts_gene
    val min_counts_cell
    val max_mito_percentage
    val symbol_col

    output:
    tuple val(meta), path("*.h5ad"), optional: true, emit: h5ad
    tuple val(meta), path("*.zarr"), optional: true, emit: zarr
    path "versions.yml"            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix     = task.ext.prefix ?: "${meta.id}_filtered"
    output_file = anndata.name.endsWith(".zarr") ? "${prefix}.zarr" : "${prefix}.h5ad"
    if (output_file == anndata.name) {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template('filter.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_filtered"
    output_file = anndata.name.endsWith(".zarr") ? "${prefix}.zarr" : "${prefix}.h5ad"
    if (output_file == anndata.name) {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    # These are needed to prevent errors during import of scanpy
    # when using singularity/apptainer
    export MPLCONFIGDIR=./tmp/mpl
    export NUMBA_CACHE_DIR=./tmp/numba

    if [[ "${output_file}" == *.zarr ]]; then
        mkdir -p "${output_file}"
    else
        touch "${output_file}"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        scanpy: \$(python3 -c 'import scanpy; print(scanpy.__version__)')
    END_VERSIONS
    """
}
