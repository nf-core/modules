process PHARMCAT_PHARMCAT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2b/2b27c134f2226e65c3be9687fdcd6dfb5eebb7998bf1ad89ff396c914fe6d81a/data' :
        'community.wave.seqera.io/library/pharmcat3:3.1.1--876b7152770ba008' }"

    input:
    tuple val(meta), path(bgz_file)

    output:
    tuple val(meta), path("*.html"),        emit: html
    tuple val(meta), path("*.json"),        emit: json
    tuple val(meta), path("*.report.json"), emit: json_input_reporter_parser
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pharmcat \\
        --matcher-vcf ${bgz_file} \\
        --base-filename ${prefix} \\
        --reporter-save-json \\
        --reporter-save-html \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: \$(pharmcat --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html
    touch ${prefix}.json
    touch ${prefix}.report.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: \$(pharmcat --version)
    END_VERSIONS
    """
}
