process PROPR_GREA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b6/b65f7192866fbd9a947df15b104808abb720e7a224bbe3ca8f7f8f680f52c97a/data' :
        'community.wave.seqera.io/library/bioconductor-limma_r-propr:f52f1d4fea746393' }"

    input:
    tuple val(meta), path(adj)
    tuple val(meta2), path(gmt)

    output:
    tuple val(meta), path("*.go.tsv"),  emit: enrichedGO
    path "versions.yml",                emit: versions
    path "*.R_sessionInfo.log",         emit: session_info

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'grea.R'
}
