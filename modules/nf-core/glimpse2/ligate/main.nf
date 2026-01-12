process GLIMPSE2_LIGATE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/glimpse-bio:2.0.1--h46b9e50_1'
        : 'biocontainers/glimpse-bio:2.0.1--h46b9e50_1'}"

    input:
    tuple val(meta), path(input_list), path(input_index)

    output:
    tuple val(meta), path("*.{vcf,bcf,vcf.gz,bcf.gz}"), emit: merged_variants
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "vcf.gz"
    """
    printf "%s\\n" ${input_list} | tr -d '[],' | sort -V > all_files.txt

    GLIMPSE2_ligate \\
        ${args} \\
        --input all_files.txt \\
        --thread ${task.cpus} \\
        --output ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glimpse2: "\$(GLIMPSE2_ligate --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]' | head -1)"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "vcf.gz"
    def create_cmd = suffix.endsWith(".gz") ? "echo | gzip > ${prefix}.${suffix}" : "touch ${prefix}.${suffix}"
    """
    ${create_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glimpse2: "\$(GLIMPSE2_ligate --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]' | head -1)"
    END_VERSIONS
    """
}
