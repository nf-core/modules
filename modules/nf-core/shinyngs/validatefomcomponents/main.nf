process SHINYNGS_VALIDATEFOMCOMPONENTS {
    tag "$sample"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-shinyngs:2.0.0--r43hdfd78af_0' :
        'biocontainers/r-shinyngs:2.0.0--r43hdfd78af_0' }"

    input:
    tuple val(meta),  path(sample), path(assay_files)
    tuple val(meta2), path(feature_meta)
    tuple val(meta3), path(contrasts)

    output:
    tuple val(meta), path("*/*.sample_metadata.tsv")   , emit: sample_meta
    tuple val(meta), path("*/*.feature_metadata.tsv")  , emit: feature_meta, optional: true
    tuple val(meta), path("*/*.assay.tsv")             , emit: assays
    tuple val(meta), path("*/*.contrasts_file.tsv")    , emit: contrasts
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // For full list of available args see
    // https://github.com/pinin4fjords/shinyngs/blob/develop/exec/validate_fom_components.R
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    def feature = feature_meta ? "--feature_metadata '$feature_meta'" : ''

    """
    validate_fom_components.R \\
        --sample_metadata "$sample" \\
        $feature \\
        --assay_files "${assay_files.join(',')}" \\
        --contrasts_file "$contrasts" \\
        --output_directory "$prefix" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-shinyngs: \$(Rscript -e "library(shinyngs); cat(as.character(packageVersion('shinyngs')))")
    END_VERSIONS
    """
}
