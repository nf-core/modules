def VERSION = '2.4.0' // PureCN outputs to stderr instead of stdout, and exits with 1 with --version

process PURECN_NORMALDB {
    tag "$meta.id"
    label 'process_medium'

    // TODO: This needs a proper container
    // cf: https://github.com/bioconda/bioconda-recipes/pull/40076
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-582ac26068889091d5e798347c637f8208d77a71:a29c64a63498b1ee8b192521fdf6ed3c65506994-0':
        'quay.io/biocontainers/mulled-v2-582ac26068889091d5e798347c637f8208d77a71:a29c64a63498b1ee8b192521fdf6ed3c65506994-0' }"


    input:
        tuple val(meta), path(coverage_files)

        val   genome
        val   assay

    output:

        tuple val(meta), path("normalDB*.rds")               , emit: rds
        tuple val(meta), path("mapping_bias*.rds")           , emit: bias_rds
        tuple val(meta), path("interval_weights*.png")       , emit: png
        tuple val(meta), path("low_coverage_targets*.png")   , emit: low_cov_png, optional: true
        path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    library_path=\$(Rscript -e 'cat(.libPaths(), sep = "\n")')
    Rscript "\$library_path"/PureCN/extdata/NormalDB.R --out-dir ./ \\
        --coverage-files $coverage_files \\
        --genome $genome \\
        --assay ${assay} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: \$(Rscript /usr/local/lib/R/library/PureCN/extdata/PureCN.R --version)
    END_VERSIONS
    """

    stub:

    """
    touch normalDB_${prefix}_${genome}.rds
    touch mapping_bias_${prefix}_${genome}.rds
    touch interval_weights_${prefix}_${genome}.png
    touch low_coverage_targets_${prefix}_${genome}.png
    """
}
