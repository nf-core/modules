process PURECN_RUN {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::bioconductor-purecn=2.4.0 bioconda::bioconductor-txdb.hsapiens.ucsc.hg38.knowngene=3.16.0 bioconductor-txdb.hsapiens.ucsc.hg19.knowngene=3.2.2 bioconda::bioconductor-org.hs.eg.db=3.16.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-582ac26068889091d5e798347c637f8208d77a71:a29c64a63498b1ee8b192521fdf6ed3c65506994-0':
        'quay.io/biocontainers/mulled-v2-582ac26068889091d5e798347c637f8208d77a71:a29c64a63498b1ee8b192521fdf6ed3c65506994-0' }"

    input:
    tuple val(meta), path(vcf)
    path tumor_coverage
    path normal_coverage
    path normal_db
    path intervals
    path mapping_bias
    path blacklist
    val genome
    val segmentation_function
    val model_type
    val seed_val

    output:
    tuple val(meta), path("*.csv")                             , emit: csv
    tuple val(meta), path("*_variants.csv")                    , emit: variants_csv
    tuple val(meta), path("*.pdf")                             , emit: pdf
    tuple val(meta), path("*.rds")                             , emit: rds
    tuple val(meta), path("*_amplification_pvalues.csv")       , emit: csv
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
    tuple val(meta), path("*-purecn-lohsummary.yaml")          , emit: loh_summary_yaml
    tuple val(meta), path("*_loh-effects.csv")                 , emit: loh_effects_csv
    tuple val(meta), path("*_segmentation.pdf")                , emit: segmentation_pdf
    tuple val(meta), path("*-sort_coverage_loess.png")         , emit: sort_coverage_loess_png
    tuple val(meta), path("*-sort_coverage_loess_qc.txt")      , emit: sort_coverage_loess_qc_txt
    tuple val(meta), path("*-sort_coverage_loess.txt.gz")      , emit: sort_coverage_loess_txt_gz
    tuple val(meta), path("*-sort_coverage.txt.gz")            , emit: sort_coverage_txt_gz
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def normal_cov = normal_coverage ? "--fun-segmentation ${normal_coverage}" : ""
    def segmentation = segmentation_function ? "--fun-segmentation ${segmentation_function}" : ""
    def map_bias = mapping_bias ? "--mapping-bias-file ${mapping_bias}" : ""
    def model = model_type ? "--model ${model_type}" : ""
    def blacklist = snp_blacklist ? "--snp-blacklist ${snp_blacklist}" : ""
    def seed = seed_val ? "--seed ${seed_val}" : ""
    def VERSION = '2.4.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    library_path=\$(Rscript -e 'cat(.libPaths(), sep = "\\n")')
    Rscript "\$library_path"PureCN.R \\
        --out ./ \\
        --tumor ${tumor_coverage} \\
        --normal ${normal_cov} \\
        --sampleid ${prefix} \\
        --vcf ${vcf} \\
        --normaldb ${normal_db} \\
        --intervals ${intervals} \\
        --genome ${genome} \\
        --stats-file ${prefix}_stats.txt \\
        ${segmentation} \\
        ${map_bias} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: ${VERSION}
    END_VERSIONS
    """

    stub:

    //TODO add STUB

    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: ${VERSION}
    END_VERSIONS
    """
}
