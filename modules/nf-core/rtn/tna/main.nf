process RTN_TNA {
    debug true
    tag "{$expression_matrix.name}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
	'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/96979cd0715edeb5d68ebbd19a353760298ebf53bf729d3b68764b2bb00683f7/data':
        'community.wave.seqera.io/library/bioconductor-rtn:2.30.0--71b797cd8b2d56b3' }"

    input:
    tuple val(meta), path(tna_object)
    tuple val(meta), path(degs)
    tuple val(meta), path(degs_log2)
    tuple val(meta), path(degs_annotation)

    output:
    tuple val(meta), path("tna.rds")               , emit: tna_object
    tuple val(meta), path("mra.rds")               , emit: mra_object
    tuple val(meta), path("gsea1.rds")             , emit: gsea1_object
    tuple val(meta), path("gsea1_plot.pdf")        , emit: gsea1_plot
    tuple val(meta), path("gsea2.rds")             , emit: gsea2_object
    tuple val(meta), path("gsea2_plot.pdf")        , emit: gsea2_plot
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    template 'rtn_tna.r'

    stub:
    def args = task.ext.args ?: ''

    """
    touch tna.rds
    touch mra.rds
    touch gsea1.rds
    touch gsea1_plot.pdf
    touch gsea2.rds
    touch gsea2_plot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-rtn: \$(Rscript -e "suppressWarnings(library(RTN)); cat(as.character(packageVersion('RTN')))")
    END_VERSIONS
    """
}
