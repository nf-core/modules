process EDGER_DIFFERENTIAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6a/6a2f2fe71499172b13396822802714ac6111e212dbfa19f189f755f52604602d/data' :
        'community.wave.seqera.io/library/bioconductor-edger:4.10.5--b2cb5dd1fa5dfb06' }"

    input:
    tuple val(meta), val(contrast_variable), val(reference), val(target), val(formula), val(comparison)
    tuple val(meta2), path(samplesheet), path(counts)

    output:
    tuple val(meta), path("*.edger.results.tsv")        , emit: results
    tuple val(meta), path("*.edger.mean_difference.png"), emit: md_plot
    tuple val(meta), path("*.edger.bcv.png")            , emit: bcv_plot
    tuple val(meta), path("*.DGEGLM.edger.rds")         , emit: rdata
    tuple val(meta), path("*.edger.model.txt")          , emit: model
    tuple val(meta), path("*.normalised_counts.tsv")    , emit: normalised_counts
    tuple val(meta), path("*.R_sessionInfo.log")        , emit: session_info
    path "versions.yml"                                 , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'edger_differential.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.edger.results.tsv
    touch ${prefix}.edger.mean_difference.png
    touch ${prefix}.edger.bcv.png
    touch ${prefix}.DGEGLM.edger.rds
    touch ${prefix}.edger.model.txt
    touch ${prefix}.normalised_counts.tsv
    touch ${prefix}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-edger: \$(Rscript -e "cat(as.character(packageVersion('edgeR')))")
    END_VERSIONS
    """
}
