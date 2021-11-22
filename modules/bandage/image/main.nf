process BANDAGE_IMAGE {
    tag "${meta.id}"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::bandage=0.8.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bandage:0.8.1--hc9558a2_2' :
        'quay.io/biocontainers/bandage:0.8.1--hc9558a2_2' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path('*.png'), emit: png
    tuple val(meta), path('*.svg'), emit: svg
    path  "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    Bandage image $gfa ${prefix}.png $args
    Bandage image $gfa ${prefix}.svg $args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(Bandage --version 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
