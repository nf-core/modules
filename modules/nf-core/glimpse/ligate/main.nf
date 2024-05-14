process GLIMPSE_LIGATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/glimpse-bio:1.1.1--hce55b13_1':
        'biocontainers/glimpse-bio:1.1.1--hce55b13_1' }"

    input:
    tuple val(meta), path(input_list), path(input_index)

    output:
    tuple val(meta), path("*.{vcf,bcf,vcf.gz,bcf.gz}"), emit: merged_variants
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "vcf.gz"
    """
    printf "%s\\n" $input_list | tr -d '[],' > all_files.txt

    GLIMPSE_ligate \\
        $args \\
        --input all_files.txt \\
        --thread $task.cpus \\
        --output ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            glimpse: "\$(GLIMPSE_ligate --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]')"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "vcf.gz"
    def args    = task.ext.args   ?: ""
    """
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glimpse: "\$(GLIMPSE_ligate --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]')"
    END_VERSIONS
    """
}
