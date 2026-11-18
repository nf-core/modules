// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CUSTOM_MTX2H5AD {

    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d1/d1a483ca814756712f25aebb7b21127890e620de38480e6d6e9e284f2d2c2e58/data'
        : 'community.wave.seqera.io/library/anndata_h5py_hdf5plugin_leidenalg_pruned:03c29b706f98b938'}"

    input:
    tuple val(meta), path(inputs), val(input_type)
    path txp2gene
    val aligner
    val kb_workflow

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    aligner = (aligner in ['cellranger', 'cellrangerarc', 'cellrangermulti']) ? 'cellranger' : aligner

    template("mtx2h5ad.py")

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.h5ad
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
