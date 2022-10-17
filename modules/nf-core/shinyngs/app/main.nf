
process SHINYNGS_APP {
    tag '$sample_sheet'
    label 'process_single'
    secret 'SHINYAPPS_TOKEN'
    secret 'SHINYAPPS_SECRET'

    conda (params.enable_conda ? "bioconda::r-shinyngs=1.2.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-shinyngs%3A1.1.0--r41hdfd78af_0':
        'quay.io/biocontainers/r-shinyngs:1.2.0--r41hdfd78af_0' }"

    input:
    path sample
    path feature_meta
    path (assay_files)
    path contrasts
    path (differential_results)

    output:
    path "app/data.rds"             , emit: data
    path "app/app.R"                , emit: app
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // For full list of available args see
    // https://github.com/pinin4fjords/shinyngs/blob/develop/exec/make_app_from_files.R
    def args = task.ext.args ?: ''

    """
    make_app_from_files.R \\
        --sample_metadata $sample \\
        --feature_metadata $feature_meta \\
        --assay_files ${assay_files.join(',')} \\
        --contrast_file $contrasts \\
        --differential_results ${differential_results.join(',')} \\
        --output_dir app \\
        --contrast_stats_assay 1 \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-shinyngs: \$(Rscript -e "library(shinyngs); cat(as.character(packageVersion('shinyngs')))")
    END_VERSIONS
    """
}
