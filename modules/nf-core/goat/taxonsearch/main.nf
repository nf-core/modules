process GOAT_TAXONSEARCH {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::goat=0.2.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/goat:0.2.0--h92d785c_0':
        'quay.io/biocontainers/goat:0.2.0--h92d785c_0' }"

    input:
    tuple val(meta), val(taxon), path(taxa_file)

    output:
    tuple val(meta), path("*.tsv"), emit: taxonsearch
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    input = taxa_file ? "-f ${taxa_file}" : "-t \"${taxon}\""
    if (!taxon && !taxa_file) error "No input. Valid input: single taxon identifier or a .txt file with identifiers"
    if (taxon && taxa_file ) error "Only one input is required: a single taxon identifier or a .txt file with identifiers"
    """
    goat-cli taxon search \\
        $args \\
        $input > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        goat: \$(goat-cli --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
