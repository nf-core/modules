process LDAK_CALCGENOTYPEERRORT2 {
    tag "${meta.id}"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
        : 'quay.io/biocontainers/r-base:4.3.1'}"

    input:
    tuple val(meta), path(he_overall_file), path(he_within_file), path(he_across_file)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: genotype_error_results
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('calc_genotype_error_t2.R')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
cat <<'TXT' > ${prefix}.txt
LDAK Genotype Error Analysis Results (T2 Statistic)
=====================================================

Stub run
TXT
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | sed -n '1s/.*version //;s/ .*//p')
    END_VERSIONS
    """
}
