process PORECHOP {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::porechop=0.2.4" : null)
    def container_image = "porechop:0.2.4--py39h7cff6ad_2"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    porechop \\
        -i $reads \\
        -t $task.cpus \\
        $args \\
        -o ${prefix}.fastq.gz \\
        > ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        porechop: \$( porechop --version )
    END_VERSIONS
    """
}
