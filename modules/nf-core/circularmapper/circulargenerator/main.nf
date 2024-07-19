// This module does the following:
//creating a modified reference genome, with an elongation_factoration of the an specified amount of bases
process CIRCULARMAPPER_CIRCULARGENERATOR {

    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/circularmapper:1.93.5--h2a3209d_3':
        'biocontainers/circularmapper:1.93.5--h2a3209d_3' }"

    input:
    tuple val(meta), path(reference)
    tuple val(meta2), val(elongation_factor)
    tuple val(meta3), val(target)

    output:
    tuple val(meta), path("*_${elongation_factor}.fasta"), emit: fasta
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    circulargenerator \
        -e ${elongation_factor} \
        -i ${reference} \
        -s ${target} \
        $args

    ## circulargenerator has a hardcoded output name. Rename if necessary to use prefix.
    if [[ "${reference.getBaseName()}_${elongation_factor}.fasta" != "${prefix}_${elongation_factor}.fasta" ]]; then
        mv ${reference.getBaseName()}_${elongation_factor}.fasta ${prefix}_${elongation_factor}.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circulargenerator: \$(circulargenerator -h | grep 'usage' | sed 's/usage: CircularGenerator//')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${elongation_factor}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circulargenerator: \$(circulargenerator -h | grep 'usage' | sed 's/usage: CircularGenerator//')
    END_VERSIONS
    """
}
