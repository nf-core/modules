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
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    quarter_reml_files_r = quarter_reml_files.sort { a, b -> a.name <=> b.name }.collect { quarterFile -> "\"${quarterFile}\"" }.join(', ')
    output_file = "${meta.id}.txt"
    template('calc_inflation.R')

    stub:
    """
    cat <<'TXT' > ${meta.id}.txt
LDAK Inflation Analysis Results
================================

Stub run
TXT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
    END_VERSIONS
    """
}
