process GOAT_TAXONSEARCH {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::goat=0.2.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/goat:0.2.0--h92d785c_0':
        'quay.io/biocontainers/goat:0.2.0--h92d785c_0' }"

    input:
    tuple val(meta), val(taxon)

    output:
    tuple val(meta), path("*.tsv"), emit: taxonsearch_results
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    goat-cli taxon search \\
        $args \\
        -t ${taxon}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        goat: \$(echo \$( goat-cli --version | cut -d' ' -f2 | head -1 2>&1) | sed 's/^.*goat //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
