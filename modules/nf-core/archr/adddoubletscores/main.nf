process ARCHR_ADDDOUBLETSCORES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "archr" // TODO update

    input:
    tuple val(meta), path(arrowfile)

    output:
    // TODO output files
    tuple val(meta), path("*.arrow")                , emit: arrowfile
    tuple val(meta), path()
    tuple val(meta), path("*.R_sessionInfo.log")    , emit: session_info
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    template adddoubletscores.r
    """

    stub:
    """
    touch QualityControl/${meta.id}/${meta.id}xx
    touch QualityControl/${meta.id}/${meta.id}xx
    touch QualityControl/${meta.id}/${meta.id}xx.pdf
    touch ${meta.id}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ArchR: \$(Rscript -e "library(ArchR); cat(as.character(packageVersion('ArchR')))")
    END_VERSIONS
    """
}
