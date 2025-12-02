process DEACON_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deacon:0.12.0--h4349ce8_0':
        'biocontainers/deacon:0.12.0--h4349ce8_0' }"

    input:
    tuple val(meta), path(index), path(reads)

    output:
    tuple val(meta), path("${prefix}*.fq.gz"), emit: fastq_filtered
    tuple val(meta), path("${prefix}.json")  , emit: log
    path "versions.yml"			             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def read_type = (reads instanceof List) ? "-o ${prefix}_1.fq -O ${prefix}_2.fq" : "> ${prefix}.fq" // deacon's automatic compression does not work
    """
    deacon \\
        filter \\
        --threads ${task.cpus} \\
        $args \\
	    --summary ${prefix}.json \\
        -d $index \\
        $reads \\
        ${read_type}

    gzip -f ${prefix}*.fq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deacon: \$(deacon --version | head -n1 | sed 's/deacon //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > '${prefix}.fq.gz'
    touch ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deacon: \$(deacon --version | head -n1 | sed 's/deacon //g')
    END_VERSIONS
    """
}
