process PURECN_RUN {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/17/171d9cdb3db28ca8a63d87dd514a97e92af353f35b8f2173173a3dc3bb801516/data':
        'community.wave.seqera.io/library/bioconductor-dnacopy_bioconductor-org.hs.eg.db_bioconductor-purecn_bioconductor-txdb.hsapiens.ucsc.hg19.knowngene_pruned:cc846801cfba58d6' }"

    input:
    tuple val(meta), path(intervals), path(coverage), path(vcf)
    path normal_db
    path mapping_bias
    val genome

    output:
    tuple val(meta), path("*.pdf")                             , emit: pdf
    tuple val(meta), path("*_local_optima.pdf")                , emit: local_optima_pdf
    tuple val(meta), path("*_dnacopy.seg")                     , emit: seg
    tuple val(meta), path("*_genes.csv")                       , emit: genes_csv                   , optional: true
    tuple val(meta), path("*_amplification_pvalues.csv")       , emit: amplification_pvalues_csv   , optional: true
    tuple val(meta), path("*.vcf.gz")                          , emit: vcf_gz                      , optional: true
    tuple val(meta), path("*_variants.csv")                    , emit: variants_csv                , optional: true
    tuple val(meta), path("*_loh.csv")                         , emit: loh_csv                     , optional: true
    tuple val(meta), path("*_chromosomes.pdf")                 , emit: chr_pdf                     , optional: true
    tuple val(meta), path("*_segmentation.pdf")                , emit: segmentation_pdf            , optional: true
    tuple val(meta), path("*_multisample.seg")                 , emit: multisample_seg             , optional: true
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_opt = vcf ? "--vcf ${vcf}": ''
    def mapping_bias_opt = mapping_bias ? "--mapping-bias-file ${mapping_bias}": ''
    def normaldb_opt = normal_db ? "--normaldb ${normal_db}": ''
    def VERSION = '2.12.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

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
        ${normaldb_opt} \\
        ${mapping_bias_opt} \\
        ${vcf_opt} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.12.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.pdf
    touch ${prefix}_local_optima.pdf
    touch ${prefix}_dnacopy.seg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: ${VERSION}
    END_VERSIONS
    """
}
