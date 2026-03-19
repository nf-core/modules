process CALC_INFLATION {
    tag "${meta.id}"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.3.1' :
        'biocontainers/r-base:4.3.1' }"

    input:
    tuple val(meta), path(ldak_reml_file), path(quarter_reml_files)

    output:
    tuple val(meta), path("${meta.id}.txt"), emit: inflation_results
    tuple val("${task.process}"), val("r-base"), eval("Rscript --version 2>&1 | cut -d\" \" -f3"), emit: versions_r_base, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def calc_inflation_script_b64 = new File("${moduleDir}/templates/calc_inflation.R").bytes.encodeBase64().toString()
    def quarter_reml_args = quarter_reml_files.sort { a, b -> a.name <=> b.name }.collect { quarterFile -> "\"${quarterFile}\"" }.join(' ')

    """
    set -euo pipefail

    printf '%s' '${calc_inflation_script_b64}' | base64 -d > calc_inflation.R

    Rscript calc_inflation.R \\
        "${ldak_reml_file}" \\
        ${quarter_reml_args} \\
        "${meta.id}.txt"
    """

    stub:
    """
    cat <<'TXT' > ${meta.id}.txt
LDAK Inflation Analysis Results
================================

Stub run
TXT
    """
}
