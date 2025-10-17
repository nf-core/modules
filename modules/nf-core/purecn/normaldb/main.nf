process PURECN_NORMALDB {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bb/bbc033d8d6415ce4883f464e2ae565df077f564da10858e4db29f6c89ad10d4a/data':
        'community.wave.seqera.io/library/bioconductor-dnacopy_bioconductor-org.hs.eg.db_bioconductor-purecn_bioconductor-txdb.hsapiens.ucsc.hg19.knowngene_pruned:7fef74d5cbdeecbe' }"


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
    def args            = task.ext.args     ?: ''
    def prefix          = task.ext.prefix   ?: "${meta.id}"
    def normal_panel    = normal_vcf        ? "--normal-panel ${normal_vcf}" : ""
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
        purecn: \$(Rscript -e 'packageVersion("PureCN")' | sed -n 's|\\[1\\] ‘\\(.*\\)’|\\1|p')
    END_VERSIONS
    """

    stub:

    def args                    = task.ext.args                     ?: ''
    def prefix                  = task.ext.prefix                   ?: "${meta.id}"
    def mapping_bias            = args.contains("--normal-panel")   ? "" : "touch mapping_bias_${prefix}_${genome}.rds"
    def mapping_bias_hq_sites   = args.contains("--normal-panel")   ? "" : "touch mapping_bias_hq_sites_${prefix}_${genome}.bed"
    """
    touch normalDB_${prefix}_${genome}.rds
    ${mapping_bias}
    ${mapping_bias_hq_sites}
    touch interval_weights_${prefix}_${genome}.png
    touch low_coverage_targets_${prefix}_${genome}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: \$(Rscript -e 'packageVersion("PureCN")' | sed -n 's|\\[1\\] ‘\\(.*\\)’|\\1|p')
    END_VERSIONS
    """
}
