process SPARSE_SIGNATURES {
    tag "$meta.id"
    label "process_long"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d410175e2fbd9c47aa685bb5dfb87cfad76d408b:a995bb98b7122825523ffed7ae131cb006e56cbe-0':
        'biocontainers/mulled-v2-d410175e2fbd9c47aa685bb5dfb87cfad76d408b:a995bb98b7122825523ffed7ae131cb006e56cbe-0' }"

    input:
    tuple val(meta), path(tsv_join,  stageAs: '*.tsv')
    val(genome)   // genome version

    output:
    tuple val(meta), path("*_mut_counts.rds"),            emit: signatures_mutCounts_rds
    tuple val(meta), path("*_cv_means_mse.rds"),          emit: signatures_cv_rds
    tuple val(meta), path("*_best_params_config.rds"),    emit: signatures_bestConf_rds
    tuple val(meta), path("*_nmf_Lasso_out.rds"),         emit: signatures_nmfOut_rds
    tuple val(meta), path("*_plot_all.rds"),              emit: signatures_plot_rds
    tuple val(meta), path("*_plot_all.pdf"),              emit: signatures_plot_pdf
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mut_counts.rds
    touch ${prefix}_cv_means_mse.rds
    touch ${prefix}_best_params_config.rds
    touch ${prefix}_nmf_Lasso_out.rds
    touch ${prefix}_plot_all.rds
    touch ${prefix}_plot_all.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-sparsesignatures: \$(Rscript -e "library(SparseSignatures); cat(as.character(packageVersion('SparseSignatures')))")
        bioconductor-bsgenome.hsapiens.1000genomes.hs37d5: \$(Rscript -e "library(BSgenome.Hsapiens.1000genomes.hs37d5); cat(as.character(packageVersion('BSgenome.Hsapiens.1000genomes.hs37d5')))")
        bioconductor-bsgenome.hsapiens.ucsc.hg38: \$(Rscript -e "library(BSgenome.Hsapiens.UCSC.hg38); cat(as.character(packageVersion('BSgenome.Hsapiens.UCSC.hg38')))")
    END_VERSIONS
    """
}
