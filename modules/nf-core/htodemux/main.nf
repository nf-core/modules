process HTODEMUX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f96b7927142847485eff858170a4cfd2d3924fb4f09de7043dd6677ac6acd09e/data':
        'community.wave.seqera.io/library/r-seurat_r-seuratobject:b11306d1bdc82827' }"

    input:
    tuple val(meta), path(seurat_object), val(assay)

    output:
    tuple val(meta), path("*_params_htodemux.csv")        , emit: params
    tuple val(meta), path("*_assignment_htodemux.csv")    , emit: assignment
    tuple val(meta), path("*_classification_htodemux.csv"), emit: classification
    tuple val(meta), path("*_htodemux.rds")               , emit: rds
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template('HTODemux.R')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_params_htodemux.csv
    touch ${prefix}_assignment_htodemux.csv
    touch ${prefix}_classification_htodemux.csv
    touch ${prefix}_htodemux.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htodemux: \$(htodemux --version)
    END_VERSIONS
    """
}
