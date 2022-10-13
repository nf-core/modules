process PYCOQC {
    tag "$summary"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pycoqc=2.5.2" : null)
        'https://depot.galaxyproject.org/singularity/pycoqc:2.5.2--py_0' :
        "${params.docker_registry ?: 'quay.io/biocontainers'}/pycoqc:2.5.2--py_0" }

    input:
    path summary

    output:
    path "*.html"        , emit: html
    path "*.json"        , emit: json
    path  "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    pycoQC \\
        $args \\
        -f $summary \\
        -o pycoqc.html \\
        -j pycoqc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pycoqc: \$(pycoQC --version 2>&1 | sed 's/^.*pycoQC v//; s/ .*\$//')
    END_VERSIONS
    """
}
