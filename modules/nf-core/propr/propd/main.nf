process PROPR_PROPD {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/27/276e08c149faa028a34ec65f74ac9eac58c2036a30db1f72acc9e35967936620/data' :
        'community.wave.seqera.io/library/bioconductor-limma_r-propr:8d299a2da993327c' }"

    input:
    tuple val(meta), val(contrast_variable), val(reference), val(target)
    tuple val(meta2), path(samplesheet), path(counts)

    output:
    tuple val(meta), path("*.propd.genewise.tsv")         , emit: results_genewise
    tuple val(meta), path("*.propd.genewise.png")         , emit: genewise_plot
    tuple val(meta), path("*.propd.rds")                  , emit: rdata                    , optional:true
    tuple val(meta), path("*.propd.pairwise.tsv")         , emit: results_pairwise         , optional:true
    tuple val(meta), path("*.propd.pairwise_filtered.tsv"), emit: results_pairwise_filtered, optional:true
    tuple val(meta), path("*.propd.adjacency.csv")        , emit: adjacency                , optional:true
    tuple val(meta), path("*.propd.fdr.tsv")              , emit: fdr                      , optional:true
    path "*.R_sessionInfo.log"                            , emit: session_info
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'propd.R'
}
