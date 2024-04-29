process GLIMPSE_CHUNK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/glimpse-bio:1.1.1--h2ce4488_2':
        'biocontainers/glimpse-bio:1.1.1--hce55b13_1' }"

    input:
    tuple val(meta), path(input), path(input_index), val(region)

    output:
    tuple val(meta), path("*.txt"), emit: chunk_chr
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def args    = task.ext.args   ?: ""

    """
    GLIMPSE_chunk \\
        $args \\
        --input $input \\
        --region $region \\
        --thread $task.cpus \\
        --output ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glimpse: "\$(GLIMPSE_chunk --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]')"
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def args    = task.ext.args   ?: ""
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glimpse: "\$(GLIMPSE_chunk --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]')"
    END_VERSIONS
    """
}
