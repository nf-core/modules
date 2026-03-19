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
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'calc_genotype_error_t2.R'

    stub:
    """
    cat <<'TXT' > ${meta.id}.txt
LDAK Genotype Error Analysis Results (T2 Statistic)
=====================================================

Stub run
TXT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
    END_VERSIONS
    """
}
