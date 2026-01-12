process PROPR_LOGRATIO {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-propr:4.2.6':
        'biocontainers/r-propr:4.2.6' }"

    input:
    tuple val(meta), path(count)

    output:
    tuple val(meta), path("*.logratio.tsv")     , emit: logratio
    tuple val(meta), path("*.R_sessionInfo.log"), emit: session_info
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'logratio.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.logratio.tsv
    touch ${prefix}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-propr: \$(Rscript -e "library('propr'); cat(as.character(packageVersion('propr')))")
    END_VERSIONS
    """
}
