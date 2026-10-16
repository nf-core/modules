process EDGER_DIFFERENTIAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/50/50c7985e75608e4fde881490b5c676741486a9b8a21fd937bcd66a39c664b47b/data'
:         'community.wave.seqera.io/library/bioconductor-edger_bioconductor-limma:9157ebf19fbdc614' }"

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
        r-base: \$(Rscript -e "cat(strsplit(version[['version.string']], ' ')[[1]][3])")
        bioconductor-edger: \$(Rscript -e "cat(as.character(packageVersion('edgeR')))")
        bioconductor-limma: \$(Rscript -e "cat(as.character(packageVersion('limma')))")
    END_VERSIONS
    """
}
