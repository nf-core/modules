process RTN_TNI {
    debug true
    tag "{$expression_matrix.name}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-rtn:2.30.0--3616043826a35de9':
        'community.wave.seqera.io/library/bioconductor-rtn:2.30.0--71b797cd8b2d56b3' }"

    input:
    tuple val(meta), path(expression_matrix)

    output:
    tuple val(meta), path("tni.rds")               , emit: tni
    tuple val(meta), path("tni_permutated.rds")    , emit: tni_perm
    tuple val(meta), path("tni_bootstrapped.rds")  , emit: tni_bootstrap
    tuple val(meta), path("tni_filtered.rds")      , emit: tni_filtered
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    template 'rtn_tni.r'

    stub:
    def args = task.ext.args ?: ''

    """
    touch tni.rds
    touch tni_permutated.rds
    touch tni_bootstrapped.rds
    touch tni_filtered.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-rtn: \$(Rscript -e "suppressWarnings(library(RTN)); cat(as.character(packageVersion('RTN')))")
    END_VERSIONS
    """
}
