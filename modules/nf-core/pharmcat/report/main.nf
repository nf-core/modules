process PHARMCAT_REPORT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://pgkb/pharmcat:3.0.1' :
        'docker.io/pgkb/pharmcat:3.0.1' }"

    input:
    tuple val(meta), path(bgz_file)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.report.json"), emit: json_input_reporter_parser
    path "versions.yml", emit: versions

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
    mkdir -p report
    echo "<html><body>Stub report</body></html>" > ${prefix}.html
    echo '{}' > ${prefix}.json
    echo '{"report": "stub"}' > ${prefix}.report.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: \$(pharmcat --version)
    END_VERSIONS
    """
}

