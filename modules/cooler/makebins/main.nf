process COOLER_MAKEBINS {
    tag "${cool_bin}"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::cooler=0.8.11" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.8.11--pyh5e36f6f_1':
        'quay.io/biocontainers/cooler:0.8.11--pyh5e36f6f_1' }"

    input:
    tuple path(chromsizes), val(cool_bin)

    output:
    path ("*.bed")       , emit: bed
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cooler makebins \\
        $args \\
        ${chromsizes} \\
        ${cool_bin} > ${prefix}_${cool_bin}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(cooler --version 2>&1 | sed 's/cooler, version //')
    END_VERSIONS
    """
}
