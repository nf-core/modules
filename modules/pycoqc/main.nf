process PYCOQC {
    tag "$summary"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pycoqc=2.5.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pycoqc:2.5.2--py_0"
    } else {
        container "quay.io/biocontainers/pycoqc:2.5.2--py_0"
    }

    input:
    path summary

    output:
    path "*.html"        , emit: html
    path "*.json"        , emit: json
    path  "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    pycoQC \\
        $args \\
        -f $summary \\
        -o pycoqc.html \\
        -j pycoqc.json

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(pycoQC --version 2>&1 | sed 's/^.*pycoQC v//; s/ .*\$//')
    END_VERSIONS
    """
}
