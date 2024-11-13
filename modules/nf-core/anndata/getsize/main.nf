process ANNDATA_GETSIZE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/04/04529821c1eff131c79f1f867fd9e8465a53ea5473bc6e4ac9405d2b9965d976/data':
        'community.wave.seqera.io/library/anndata:0.10.9--1eab54e300e1e584' }"

    input:
    tuple val(meta), path(h5ad)
    val size_type

    output:
    tuple val(meta), path("*.txt"), emit: size
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'getsize.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        anndata: \$(python3 -c 'import anndata; print(anndata.__version__)')
    END_VERSIONS
    """
}
