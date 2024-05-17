// This module does the following:
//creating a modified reference genome, with an elongation of the an specified amount of bases
process CIRCULARMAPPER_CIRCULARGENERATOR {

    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/circularmapper:1.93.5--h2a3209d_3':
        'biocontainers/circularmapper:1.93.5--h2a3209d_3' }"

    input:
    tuple val(meta), path(reference)
    val(elong)

    output:
    tuple val(meta), path("*_${elong}.fasta"), emit: fasta
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    circulargenerator -e ${elong} \
        -i ${reference} \
        -s ${prefix} \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circulargenerator: \$(circulargenerator -h | grep 'usage' | sed 's/usage: CircularGenerator//')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${elong}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circulargenerator: \$(circulargenerator -h | grep 'usage' | sed 's/usage: CircularGenerator//')
    END_VERSIONS
    """
}
