process SOURMASH_TAXANNOTATE {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::sourmash=4.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.5.0--hdfd78af_0':
        'quay.io/biocontainers/sourmash:4.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(gather_results)
    path(taxonomy)

    output:
    tuple val(meta), path("*.with-lineages.csv.gz"), emit: result
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    # Currently, `sourmash tax annotate` does not support gz-compressed input,
    # streaming is also not supported
    # so we need to decomress the input file first
    gunzip -c ${gather_results} > ${prefix}.csv

    sourmash \\
        tax annotate \
        $args \\
        --gather-csv ${prefix}.csv \\
        --taxonomy ${taxonomy} \\
        --output-dir "."

    ## Output file name = ${prefix}.with-lineages.csv

    ## Compress output
    gzip "${prefix}.with-lineages.csv"

    # Remove decompressed file
    rm "${prefix}.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}
