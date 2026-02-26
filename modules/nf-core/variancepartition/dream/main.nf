process VARIANCEPARTITION_DREAM {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d6/d6fa8a7908dd357484d63bb574a378f1739440a2c0752ce91b1e0e9d1ac1638c/data' :
        'community.wave.seqera.io/library/bioconductor-edger_bioconductor-variancepartition_r-optparse:ba778938d72f30c5' }"

    input:
    tuple val(meta), val(contrast_variable), val(reference), val(target), val(formula), val(comparison)
    tuple val(meta2), path(samplesheet), path(counts)

    output:
    tuple val(meta), path("*.dream.results.tsv")        , emit: results
    tuple val(meta), path("*.dream.model.txt")          , emit: model
    tuple val(meta), path("*.normalised_counts.tsv")    , emit: normalised_counts, optional: true
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'dream.R'

    stub:
    """
    touch "${meta.id}.dream.results.tsv"
    touch "${meta.id}.dream.model.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-edger: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
        r-variancepartition: \$(Rscript -e "library(variancePartition); cat(as.character(packageVersion('variancePartition')))")
    END_VERSIONS
    """
}
