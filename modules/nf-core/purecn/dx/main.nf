process PURECN_RUN {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-582ac26068889091d5e798347c637f8208d77a71:a29c64a63498b1ee8b192521fdf6ed3c65506994-0':
        'biocontainers/mulled-v2-582ac26068889091d5e798347c637f8208d77a71:a29c64a63498b1ee8b192521fdf6ed3c65506994-0' }"

    input:
    tuple val(meta), path(rds_file)
    path callable_regions
    path exclude_regions

    output:
    tuple val(meta), path("*_mutation_burden_csv")             , emit: mutation_burden_csv
    tuple val(meta), path("*_cin_csv")                         , emit: cin_csv
    tuple val(meta), path("*_signatures_csv")                  , emit: signatures_csv              , optional: true
    tuple val(meta), path("*_signatures_pdf")                  , emit: signatures_pdf              , optional: true
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.4.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    library_path=\$(Rscript -e 'cat(.libPaths(), sep = "\\n")')
    Rscript "\$library_path"/PureCN/extdata/Dx.R \\
        --rds ${rds_file} \\
        --callable ${callable_regions} \\
        --exclude${exclude_regions} \\
        --out ./
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
    touch ${prefix}.pdf
    touch ${prefix}_local_optima.pdf
    touch ${prefix}_dnacopy.seg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: ${VERSION}
    END_VERSIONS
    """
}
