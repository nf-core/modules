process LEVIOSAM2_INDEX {
    tag "$meta.id"
    label 'process_low'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/leviosam2:0.4.2--h4ac6f70_0':
        'biocontainers/leviosam2:0.4.2--h4ac6f70_0' }"

    input:

    tuple val(meta), path(fai)
    path(chain)

    output:
    tuple val(meta), path("*.clft"), emit: clft
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    leviosam2 \\
        index \\
        -c ${chain} \\
        -p ${prefix} \\
        -F ${fai}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leviosam2: \$(leviosam2 --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.clft

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leviosam2: \$(leviosam2 --version)
    END_VERSIONS
    """
}

    