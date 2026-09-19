process PROPR_PROPR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-propr:5.0.3':
        'quay.io/biocontainers/r-propr:5.0.3' }"

    input:
    tuple val(meta), path(count)

    output:
    tuple val(meta), path("*.propr.rds"), emit: propr
    tuple val(meta), path("*.propr.tsv"), emit: matrix
    tuple val(meta), path("*.fdr.tsv"),   emit: fdr         , optional:true
    tuple val(meta), path("*.adj.csv"),   emit: adj         , optional:true
    path "*.warnings.log",                emit: warnings
    path "*.R_sessionInfo.log",           emit: session_info
    path "versions.yml",                  emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'propr.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.propr.rds
    touch ${prefix}.propr.tsv
    touch ${prefix}.fdr.tsv
    touch ${prefix}.adj.csv
    touch ${prefix}.warnings.log
    touch ${prefix}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-propr: \$(Rscript -e "cat(as.character(packageVersion('propr')))")
    END_VERSIONS
    """
}
