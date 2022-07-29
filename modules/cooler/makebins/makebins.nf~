process COOLER_MAKEBINS {
    tag "${cool_bin}"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::cooler=0.8.11" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.8.11--pyh3252c3a_0' :
        'quay.io/biocontainers/cooler:0.8.11--pyh3252c3a_0' }"

    input:
    tuple path(chromsizes), val(cool_bin)

    output:
    path ("*.bed")       , emit: bed
    path ("versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    cooler makebins \\
        $args \\
        ${chromsizes} \\
        ${cool_bin} > cooler_bins_${cool_bin}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(cooler --version 2>&1 | sed 's/cooler, version //')
    END_VERSIONS
    """
}
