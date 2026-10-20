process CUSTOM_MTX2H5AD {

    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/74/74d1b20458a4cdb5ecbe3ead3bca69950d3d5608747274c0ee02a30095ba9ac6/data':
        'community.wave.seqera.io/library/hdf5plugin_leidenalg_python-igraph_scanpy:5ac4f1cde98e0026' }"

    input:
    tuple val(meta), path(inputs), val(input_type)
    path txp2gene
    val aligner
    val kb_workflow

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "versions.yml"                    , emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    aligner = (aligner in [ 'cellranger', 'cellrangerarc', 'cellrangermulti' ]) ? 'cellranger' : aligner

    template "mtx2h5ad.py"

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
