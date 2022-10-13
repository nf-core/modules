process BEDTOOLS_SPLIT {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    def container_image = "/bedtools:2.30.0--h468198e_3"
                                                 container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(bed)
    val(number_of_files)

    output:
    tuple val(meta), path("*.bed"), emit: beds
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bedtools \\
        split \\
        $args \\
        -i $bed \\
        -p $prefix \\
        -n $number_of_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
