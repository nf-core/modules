process PURECN_NORMALDB {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::bioconductor-purecn=2.4.0 bioconda::bioconductor-txdb.hsapiens.ucsc.hg38.knowngene=3.16.0 bioconductor-txdb.hsapiens.ucsc.hg19.knowngene=3.2.2 bioconda::bioconductor-org.hs.eg.db=3.16.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-582ac26068889091d5e798347c637f8208d77a71:a29c64a63498b1ee8b192521fdf6ed3c65506994-0':
        'quay.io/biocontainers/mulled-v2-582ac26068889091d5e798347c637f8208d77a71:a29c64a63498b1ee8b192521fdf6ed3c65506994-0' }"

    input:
    tuple val(meta), path(coverage_files), path(normal_vcf), path(normal_vcf_tbi)
    val   genome
    val   assay

    output:
    tuple val(meta), path("normalDB*.rds")               , emit: rds
    tuple val(meta), path("interval_weights*.png")       , emit: png
    tuple val(meta), path("mapping_bias*.rds")           , emit: bias_rds,    optional: true
    tuple val(meta), path("mapping_bias_hq_sites*.bed")  , emit: bias_bed,    optional: true
    tuple val(meta), path("low_coverage_targets*.bed")   , emit: low_cov_bed, optional: true
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def normal_panel = normal_vcf ? "--normal-panel ${normal_vcf}" : ""
    def VERSION = '2.4.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    echo $coverage_files | tr ' ' '\\n' > coverages.list
    library_path=\$(Rscript -e 'cat(.libPaths(), sep = "\\n")')
    Rscript "\$library_path"/PureCN/extdata/NormalDB.R --out-dir ./ \\
        --coverage-files coverages.list \\
        --genome ${genome} \\
        --assay ${assay} \\
        ${normal_panel} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: ${VERSION}
    END_VERSIONS
    """

    stub:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mapping_bias = args.contains("--normal-panel") ? "" : "touch mapping_bias_${prefix}_${genome}.rds"
    def mapping_bias_hq_sites = args.contains("--normal-panel") ? "" : "touch mapping_bias_hq_sites_${prefix}_${genome}.bed"
    def VERSION = '2.4.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch normalDB_${prefix}_${genome}.rds
    ${mapping_bias}
    ${mapping_bias_hq_sites}
    touch interval_weights_${prefix}_${genome}.png
    touch low_coverage_targets_${prefix}_${genome}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: ${VERSION}
    END_VERSIONS
    """
}
