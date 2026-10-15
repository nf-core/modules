process DEXSEQ_DTU {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e4/e42feb37a257a80b9a48a628d05a6f3a1283bed077b0bb0a6a97bb7a57c25c5f/data'
        : 'community.wave.seqera.io/library/bioconductor-dexseq:1.56.0--548bf2f0349db29a'}"

    input:
    tuple val(meta), val(contrast_variable), val(reference), val(target)
    tuple val(meta2), path(samplesheet), path(counts)

    output:
    tuple val(meta), path("*.dexseq.results.tsv")        , emit: results
    tuple val(meta), path("*.dexseq.gene_qvalues.tsv")   , emit: gene_qvalues
    tuple val(meta), path("*.dexseq.dataset.rds")        , emit: dataset_rds
    tuple val(meta), path("*.dexseq.results.rds")        , emit: results_rds
    tuple val(meta), path("*.dexseq.gene_qvalues.rds")   , emit: gene_qvalues_rds
    tuple val(meta), path("*.R_sessionInfo.log")         , emit: session_info
    path "versions.yml"                                  , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'dexseq_dtu.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dexseq.results.tsv
    touch ${prefix}.dexseq.gene_qvalues.tsv
    touch ${prefix}.dexseq.dataset.rds
    touch ${prefix}.dexseq.results.rds
    touch ${prefix}.dexseq.gene_qvalues.rds
    touch ${prefix}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(paste(R.version[['major']], R.version[['minor']], sep = '.'))")
        bioconductor-dexseq: \$(Rscript -e "cat(as.character(packageVersion('DEXSeq')))")
    END_VERSIONS
    """
}
