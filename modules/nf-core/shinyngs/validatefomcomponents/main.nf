process SHINYNGS_VALIDATEFOMCOMPONENTS {
    tag '$sample'
    label 'process_single'

    conda (params.enable_conda ? "bioconda::r-shinyngs=1.4.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-shinyngs%3A1.4.1--r41hdfd78af_0':
        'quay.io/biocontainers/r-shinyngs:1.4.1--r41hdfd78af_0' }"

    input:
    tuple val(meta), path(sample), path(feature_meta), path(assay_files)
    tuple val(meta2), path(contrasts)

    output:
    tuple val(meta), path("*/*.sample_metadata.tsv"), path("*/*.feature_metadata.tsv"), path("*/*.assay.tsv")   , emit: fom
    path("*/*.contrasts_file.tsv")                                                                              , emit: contrasts
    path "versions.yml"                                                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // For full list of available args see
    // https://github.com/pinin4fjords/shinyngs/blob/develop/exec/validate_fom_components.R
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id

    """
    validate_fom_components.R \\
        --sample_metadata $sample \\
        --feature_metadata $feature_meta \\
        --assay_files ${assay_files.join(',')} \\
        --contrasts_file $contrasts \\
        --output_directory $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-shinyngs: \$(Rscript -e "library(shinyngs); cat(as.character(packageVersion('shinyngs')))")
    END_VERSIONS
    """
}
