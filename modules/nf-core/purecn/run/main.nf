process PURECN_RUN {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::bioconductor-purecn=2.4.0 bioconda::bioconductor-txdb.hsapiens.ucsc.hg38.knowngene=3.16.0 bioconductor-txdb.hsapiens.ucsc.hg19.knowngene=3.2.2 bioconda::bioconductor-org.hs.eg.db=3.16.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-582ac26068889091d5e798347c637f8208d77a71:a29c64a63498b1ee8b192521fdf6ed3c65506994-0':
        'biocontainers/mulled-v2-582ac26068889091d5e798347c637f8208d77a71:a29c64a63498b1ee8b192521fdf6ed3c65506994-0' }"

    input:
    tuple val(meta), path(intervals), path(coverage)
    path normal_db
    val genome

    output:
    tuple val(meta), path("*.csv")                             , emit: csv
    tuple val(meta), path("*_variants.csv")                    , emit: variants_csv
    tuple val(meta), path("*.pdf")                             , emit: pdf
    tuple val(meta), path("*.rds")                             , emit: rds
    tuple val(meta), path("*_amplification_pvalues.csv")       , emit: amplification_pvalues_csv
    tuple val(meta), path("*_chromosomes.pdf")                 , emit: chr_pdf
    tuple val(meta), path("*_dnacopy.seg")                     , emit: seg
    tuple val(meta), path("*_genes.csv")                       , emit: genes_csv
    tuple val(meta), path("*_local_optima.pdf")                , emit: local_optima_pdf
    tuple val(meta), path("*.log")                             , emit: log
    tuple val(meta), path("*_loh.vcf")                         , emit: loh_vcf
    tuple val(meta), path("*_loh.vcf.gz")                      , emit: loh_vcf_gz
    tuple val(meta), path("*_loh.vcf.gz.tbi")                  , emit: loh_vcf_tbi
    tuple val(meta), path("*_loh.csv")                         , emit: loh_csv
    tuple val(meta), path("*_loh-effects-stats.csv")           , emit: loh_stats_csv
    tuple val(meta), path("*_loh-effects-stats.genes.txt")     , emit: loh_effects_txt
    tuple val(meta), path("*_loh-effects-stats.html")          , emit: loh_effects_html
    tuple val(meta), path("*_loh-effects.vcf.gz")              , emit: loh_effects_vcf
    tuple val(meta), path("*_loh-effects.vcf.gz.tbi")          , emit: loh_effects_tbi
    tuple val(meta), path("*_loh-summary.yaml")                , emit: loh_summary_yaml
    tuple val(meta), path("*_loh-effects.csv")                 , emit: loh_effects_csv
    tuple val(meta), path("*_segmentation.pdf")                , emit: segmentation_pdf
    tuple val(meta), path("*_sort_coverage_loess.png")         , emit: sort_coverage_loess_png
    tuple val(meta), path("*_sort_coverage_loess_qc.txt")      , emit: sort_coverage_loess_qc_txt
    tuple val(meta), path("*_sort_coverage_loess.txt.gz")      , emit: sort_coverage_loess_txt_gz
    tuple val(meta), path("*_sort_coverage.txt.gz")            , emit: sort_coverage_txt_gz
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.4.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    library_path=\$(Rscript -e 'cat(.libPaths(), sep = "\\n")')
    Rscript "\$library_path"/PureCN/extdata/PureCN.R \\
        --out ./ \\
        --tumor ${coverage} \\
        --sampleid ${prefix} \\
        --normaldb ${normal_db} \\
        --intervals ${intervals} \\
        --genome ${genome} \\
        --parallel \\
        --cores ${task.cpus} \\
        --stats-file ${prefix}_stats.txt \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.4.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.csv
    touch ${prefix}_variants.csv
    touch ${prefix}.pdf
    touch ${prefix}.rds
    touch ${prefix}_amplification_pvalues.csv
    touch ${prefix}_chromosomes.pdf
    touch ${prefix}_dnacopy.seg
    touch ${prefix}_genes.csv
    touch ${prefix}_local_optima.pdf
    touch ${prefix}.log
    touch ${prefix}_loh.vcf
    touch ${prefix}_loh.vcf.gz
    touch ${prefix}_loh.vcf.gz.tbi
    touch ${prefix}_loh.csv
    touch ${prefix}_loh-effects-stats.csv
    touch ${prefix}_loh-effects-stats.genes.txt
    touch ${prefix}_loh-effects-stats.html
    touch ${prefix}_loh-effects.vcf.gz
    touch ${prefix}_loh-effects.vcf.gz.tbi
    touch ${prefix}_loh-summary.yaml
    touch ${prefix}_loh-effects.csv
    touch ${prefix}_segmentation.pdf
    touch ${prefix}_sort_coverage_loess.png
    touch ${prefix}_sort_coverage_loess_qc.txt
    touch ${prefix}_sort_coverage_loess.txt.gz
    touch ${prefix}_sort_coverage.txt.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: ${VERSION}
    END_VERSIONS
    """
}
