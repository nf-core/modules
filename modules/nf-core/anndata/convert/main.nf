process ANNDATA_CONVERT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d8036a262b63ab783572e3cb3120b6b138593e19d93a1059b9bfef80785c7218/data' :
        'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-rhdf5_bioconductor-singlecellexperiment_r-seurat:a0f51df063bb9b2a' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${prefix}.seurat.rds"), emit: seurat
    tuple val(meta), path("${prefix}.sce.rds"), emit: sce
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    template 'convert.R'

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.seurat.rds
    touch ${prefix}.sce.rds
    touch versions.yml
    """
}
