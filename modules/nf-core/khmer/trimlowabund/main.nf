process KHMER_TRIMLOWABUND {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/khmer:3.0.0a3--py37haa7609a_2' :
        'biocontainers/khmer:3.0.0a3--py37haa7609a_2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.abundtrim"), emit: sequence
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if (!task.memory) {
        log.info '[KHMER_TRIMLOWABUND] Available memory not known - defaulting to 16GB. Specify process memory requirements to change this.'
        avail_mem = 16
    } else {
        avail_mem = task.memory.toGiga()
    }
    """
    trim-low-abund.py \\
        -M ${avail_mem}e9 \\
        $args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        khmer: \$( trim-low-abund.py --version 2>&1 | grep ^khmer | sed 's/^khmer //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta}.abundtrim

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        khmer: \$( trim-low-abund.py --version 2>&1 | grep ^khmer | sed 's/^khmer //' )
    END_VERSIONS
    """
}
