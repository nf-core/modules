process GLIMPSE_LIGATE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::glimpse-bio=1.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/glimpse-bio:1.1.1--hce55b13_1':
        'quay.io/biocontainers/glimpse-bio:1.1.1--hce55b13_1' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.{vcf,bcf}"), emit: merged_variant
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "vcf"
    """
    GLIMPSE_ligate \\
        $args \\
        --input $input \\
        --thread $task.cpus \\
        --output ${prefix}_merged.${suffix}

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            glimpse: "\$(GLIMPSE_ligate --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]')"
    END_VERSIONS
    """
}
