process CALC_GENOTYPE_ERROR_T2 {
    tag "${meta.id}"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.3.1' :
        'biocontainers/r-base:4.3.1' }"

    input:
    tuple val(meta), path(he_overall_file), path(he_within_file), path(he_across_file)

    output:
    tuple val(meta), path("${meta.id}.txt"), emit: genotype_error_results
    tuple val("${task.process}"), val("r-base"), eval("Rscript --version 2>&1 | cut -d\" \" -f3"), emit: versions_r_base, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def calc_genotype_error_t2_script_b64 = new File("${moduleDir}/templates/calc_genotype_error_t2.R").bytes.encodeBase64().toString()

    """
    set -euo pipefail

    printf '%s' '${calc_genotype_error_t2_script_b64}' | base64 -d > calc_genotype_error_t2.R

    Rscript calc_genotype_error_t2.R \\
        "${he_overall_file}" \\
        "${he_within_file}" \\
        "${he_across_file}" \\
        "${meta.id}.txt"
    """

    stub:
    """
    cat <<'TXT' > ${meta.id}.txt
LDAK Genotype Error Analysis Results (T2 Statistic)
=====================================================

Stub run
TXT
    """
}
