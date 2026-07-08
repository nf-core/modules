process LDAK_CALCINFLATION {
    tag "${meta.id}"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
        : 'quay.io/biocontainers/r-base:4.3.1'}"

    input:
    tuple val(meta), path(ldak_reml_file), path(quarter_reml_files)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: inflation_results
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    quarter_reml_files_r = quarter_reml_files.sort { a, b -> a.name <=> b.name }.collect { quarterFile -> "\"${quarterFile}\"" }.join(', ')
    template('calc_inflation.R')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
cat <<'TXT' > ${prefix}.txt
LDAK Inflation Analysis Results
================================

Stub run
TXT
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | sed -n '1s/.*version //;s/ .*//p')
    END_VERSIONS
    """
}
