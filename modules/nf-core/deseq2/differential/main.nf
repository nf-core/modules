process DESEQ2_DIFFERENTIAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a15f5d61792b60b6179afd885db27d3fe60eb4c42e805c8887ed0416d88cb484/data' :
        'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-limma:b56a0c9ddc3e87e1' }"

    input:
    tuple val(meta), val(contrast_variable), val(reference), val(target), val(formula), val(comparison)
    tuple val(meta2), path(samplesheet), path(counts)
    tuple val(meta3), path(control_genes_file)
    tuple val(meta4), path(transcript_lengths_file)

    output:
    tuple val(meta), path("*.deseq2.results.tsv")              , emit: results
    tuple val(meta), path("*.deseq2.dispersion.png")           , emit: dispersion_plot_png
    tuple val(meta), path("*.deseq2.dispersion.pdf")           , emit: dispersion_plot_pdf
    tuple val(meta), path("*.dds.rld.rds")                     , emit: rdata
    tuple val(meta), path("*.deseq2.sizefactors.tsv")          , emit: size_factors
    tuple val(meta), path("*.normalised_counts.tsv")           , emit: normalised_counts
    tuple val(meta), path("*.rlog.tsv")                        , optional: true, emit: rlog_counts
    tuple val(meta), path("*.vst.tsv")                         , optional: true, emit: vst_counts
    tuple val(meta), path("*.deseq2.model.txt")                , emit: model
    tuple val(meta), path("*.R_sessionInfo.log")               , emit: session_info
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'deseq2_differential.R'

    stub:
    """
    touch ${meta.id}.deseq2.results.tsv
    touch ${meta.id}.deseq2.dispersion.png
    touch ${meta.id}.deseq2.dispersion.pdf
    touch ${meta.id}.dds.rld.rds
    touch ${meta.id}.deseq2.sizefactors.tsv
    touch ${meta.id}.normalised_counts.tsv
    touch ${meta.id}.rlog.tsv
    touch ${meta.id}.deseq2.model.txt
    touch ${meta.id}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}
