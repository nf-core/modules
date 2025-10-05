process MSSTATS_MSSTATSLFQ {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-msstats:4.14.0--r44he5774e6_0' :
        'biocontainers/bioconductor-msstats:4.14.0--r44he5774e6_0' }"

    input:
    tuple val(meta), path(msstats_csv_input)

    output:
    // The generation of the PDFs from MSstats are very unstable, especially with auto-contrasts.
    // And users can easily fix anything based on the csv and the included script -> make optional
    tuple val(meta), path("*.pdf"), emit: pdf, optional: true
    tuple val(meta), path("*.csv"), emit: msstats_csv
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    msstats_plfq.R \\
        ${msstats_csv_input} \\
        ${args} \\
        2>&1 | tee msstats.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-msstats: \$(Rscript -e "library(MSstats); cat(as.character(packageVersion('MSstats')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    touch test_sample_msstats_in_msstats_results.csv
    touch msstats.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-msstats: \$(Rscript -e "library(MSstats); cat(as.character(packageVersion('MSstats')))")
    END_VERSIONS
    """
}
