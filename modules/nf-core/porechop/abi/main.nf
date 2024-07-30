process PORECHOP_ABI {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/porechop_abi:0.5.0--py310h590eda1_0':
        'biocontainers/porechop_abi:0.5.0--py310h590eda1_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz") , emit: reads
    tuple val(meta), path("*.log")      , emit: log
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.porechop_abi"
    if ("$reads" == "${prefix}.fastq.gz") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    porechop_abi \\
        --input $reads \\
        --threads $task.cpus \\
        $args \\
        --output ${prefix}.fastq.gz \\
        | tee ${prefix}.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        porechop_abi: \$( porechop_abi --version )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.porechop_abi"
    """
    echo "" | gzip > ${prefix}.fastq.gz
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        porechop_abi: \$( porechop_abi --version )
    END_VERSIONS
    """
}
