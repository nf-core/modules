process SOURMASH_TAXANNOTATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.8.14--hdfd78af_0':
        'biocontainers/sourmash:4.8.14--hdfd78af_0' }"

    input:
    tuple val(meta), path(gather_results)
    path(taxonomy)

    output:
    tuple val(meta), path("*.with-lineages.csv.gz"), emit: result
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    sourmash \\
        tax annotate \
        $args \\
        --gather-csv ${gather_results} \\
        --taxonomy ${taxonomy} \\
        --output-dir "."

    ## Compress output
    gzip --no-name *.with-lineages.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}
