process ARCHR_CREATEARROWFILES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "archr_dckr" // TODO update

    input:
    tuple val(meta), path(fragments)

    output:
    tuple val(meta), path("*.arrow")                , emit: arrowfile
    tuple val(meta), path("*.rds")                  , emit: arrowmetadata
    tuple val(meta), path("*.pdf")                   , emit: qualitycontrolplots
    tuple val(meta), path("*.R_sessionInfo.log")    , emit: session_info
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    template 'createarrowfiles.r'


    stub:
    """
    touch ${meta.id}.arrow
    touch ${meta.id}-Pre-Filter-Metadata.rds
    touch ${meta.id}-Fragment_Size_Distribution.pdf
    touch ${meta.id}-TSS_by_Unique_Frags.pdf
    touch ${meta.id}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ArchR: \$(Rscript -e "library(ArchR); cat(as.character(packageVersion('ArchR')))")
    END_VERSIONS
    """
}
