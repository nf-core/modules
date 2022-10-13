process GENOTYPHI_PARSE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::genotyphi=1.9.1" : null)
    def container_image = "/genotyphi:1.9.1--hdfd78af_1"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }
n

    input:
    tuple val(meta), path(json)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_typhi_mykrobe.py \\
        --jsons $json \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genotyphi: \$(echo \$(genotyphi --version 2>&1) | sed 's/^.*GenoTyphi v//;' )
    END_VERSIONS
    """
}
