process SCANPY_HASHSOLO {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/168ecbbe27ccef766741ccbf937b0d2675be2e19b0565035e0719f1e9ea5ee95/data':
        'community.wave.seqera.io/library/python_pyyaml_scanpy:b5509a698e9aae25' }"

    input:
    tuple val(meta), path(hto_data)

    output:
    tuple val(meta), path("*_hashsolo.csv"), emit: assignment
    tuple val(meta), path("*_hashsolo.jpg"), emit: plot
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    use_rna_data                  = task.ext.use_rna_data                  ?: false    // Whether to use RNA data for demultiplexing
    priors                        = task.ext.priors_negative               ?: [0.01, 0.8, 0.19]  // A list of prior for each hypothesis. The first element is prior for the negative hypothesis, second for the singlet hypothesis, third element for the doublet hypothesis
    pre_existing_clusters         = task.ext.pre_existing_clusters         ?: null            // Column in cell_hashing_adata.obs for how to break up demultiplexing
    number_of_noise_barcodes      = task.ext.number_of_noise_barcodes      ?: null     // Number of barcodes to use to create noise distribution
    clustering_data               = task.ext.clustering_data               ?: null     // Input directory containing transcriptomic data in 10x mtx format

    template('hashsolo.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad

    # Prevent failures during scanpy import
    export MPLCONFIGDIR=./tmp/mpl
    export NUMBA_CACHE_DIR=./tmp/numba

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -f 2 -d " ")
        scanpy: \$(python3 -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
