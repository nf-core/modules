process BEDTOOLS_MAKEWINDOWS {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
        def container_image = "/bedtools:2.30.0--h7d7f7ad_1"
                                                       container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(regions)
    val(use_bed)

    output:
    tuple val(meta), path("*.tab"), emit: tab
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def arg_input = use_bed ? "-b $regions" : "-g $regions"
    """
    bedtools \\
        makewindows \\
        ${arg_input} \\
        $args \\
        > ${prefix}.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
