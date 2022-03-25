process CONTROLFREEC_MAKEGRAPH {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::control-freec=11.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/control-freec:11.6--h1b792b2_1':
        'quay.io/biocontainers/control-freec:11.6--h1b792b2_1' }"

    input:
    tuple val(meta), path(ratio), path(baf)

    output:
    tuple val(meta), path("*_BAF.txt.png")       , emit: png_baf
    tuple val(meta), path("*_ratio.txt.log2.png"), emit: png_ratio_log2
    tuple val(meta), path("*_ratio.txt.png")     , emit: png_ratio

    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def baf = baf ?: ""
    """
    cat /usr/local/bin/makeGraph.R | R --slave --args ${args} ${ratio} ${baf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: \$(echo \$(freec -version 2>&1) | sed 's/^.*Control-FREEC  //; s/:.*\$//' | sed -e "s/Control-FREEC v//g" )
    END_VERSIONS
    """
}
