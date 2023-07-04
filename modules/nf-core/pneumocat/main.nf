VERSION = '1.2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
process PNEUMOCAT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::pneumocat=1.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pneumocat:1.2.1--0':
        'biocontainers/pneumocat:1.2.1--0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.xml"), emit: xml
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    PneumoCaT.py \\
        --input_directory ./ \\
        $args \\
        --threads $task.cpus \\
        --output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pneumocat: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.xml
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pneumocat: $VERSION
    END_VERSIONS
    """
}
