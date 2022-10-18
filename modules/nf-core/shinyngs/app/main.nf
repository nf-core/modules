
process SHINYNGS_APP {
    tag '$sample_sheet'
    label 'process_single'

    // To be able to pass the necessary secrets for shinyapps.io deployment,
    // this process must be configured by placing a statement like the
    // following in the nextflow.config:
    //
    // withName: SHINYNGS_APP {
    //     secret 'SHINYAPPS_TOKEN'
    //     secret 'SHINYAPPS_SECRET
    // }
    //
    // Those values must then be set in your Nextflow secrets.

    conda (params.enable_conda ? "bioconda::r-shinyngs=1.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-shinyngs%3A1.3.0--r41hdfd78af_0':
        'quay.io/biocontainers/r-shinyngs:1.3.0--r41hdfd78af_0' }"

    input:
    tuple val(meta), path(sample)
    tuple val(meta), path(feature_meta)
    tuple val(meta), path(assay_files)
    tuple val(meta), path(contrasts)
    tuple val(meta), path(differential_results)

    output:
    tuple val(meta), path("*/data.rds"), path("*/app.R")    , emit: data
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // For full list of available args see
    // https://github.com/pinin4fjords/shinyngs/blob/develop/exec/make_app_from_files.R
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id

    """
    make_app_from_files.R \\
        --sample_metadata $sample \\
        --feature_metadata $feature_meta \\
        --assay_files ${assay_files.join(',')} \\
        --contrast_file $contrasts \\
        --differential_results ${differential_results.join(',')} \\
        --output_dir $prefix \\
        --contrast_stats_assay 1 \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-shinyngs: \$(Rscript -e "library(shinyngs); cat(as.character(packageVersion('shinyngs')))")
    END_VERSIONS
    """
}
