process PYCHOPPER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pychopper:2.7.10--pyhdfd78af_0':
        'biocontainers/pychopper:2.7.10--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.out.fastq.gz"), emit: fastq
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    pychopper \\
        $args \\
        -t $task.cpus \\
        $fastq \\
        ${prefix}.out.fastq \\

    gzip ${prefix}.out.fastq > ${prefix}.out.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pychopper: 2.7.10 (hard coded- check container used for this module)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.out.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pychopper: 2.7.10 (hard coded- check container used for this module)
    END_VERSIONS
    """
}
