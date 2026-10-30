process CUSTOM_MTX2H5AD {

    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/02/029644f69918a0d2abdc7a5438decb054e2a4b292e87536afc5589cceff581d8/data'
        : 'community.wave.seqera.io/library/hdf5plugin_leidenalg_python-igraph_python_scanpy:20f9779bbe32a4c8'}"

    input:
    tuple val(meta), path(inputs), val(input_type)
    path txp2gene
    val aligner
    val kb_workflow

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

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
