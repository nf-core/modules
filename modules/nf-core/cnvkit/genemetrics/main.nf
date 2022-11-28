
process CNVKIT_GENEMETRICS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::cnvkit=0.9.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.9--pyhdfd78af_0':
        'quay.io/biocontainers/cnvkit:0.9.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(cnr), path(cns)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    //tuple val(meta), path("*.cnn"), emit: cnn
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def segments = cns ? '--segment $cns' : ''

    """
    cnvkit.py \\
        genemetrics \\
        $cnr \\
        $segments \\
        --output ${prefix}.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed -e "s/cnvkit v//g")
    END_VERSIONS
    """
}
