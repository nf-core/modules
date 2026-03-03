process EIGSCIS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooltools:0.7.1--py39hff726c5_2':
        'biocontainers/cooltools:0.7.1--py39hff726c5_2' }"

    input:
    tuple val(meta), path(cool)
    path(bed)

    output:
    tuple val(meta), path("*compartment*"), emit:result
    tuple val(meta), path("*.bw"), emit: bigwig, optional: true
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def phasing_track = bed ? "--phasing-track ${bed}" : ""

    """
    cooltools \\
        eigs-cis ${args} \\
        $phasing_track \\
        -o ${prefix}_compartments.bed ${cool}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooltools: \$(cooltools --version | grep 'cooltools, version ' | sed 's/cooltools, version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_compartments.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooltools: \$(cooltools --version 2>&1 | grep version | sed 's/cooltools, version //')
    END_VERSIONS
    """

}
