process UPDHMM_COLLAPSEEVENTS {
    tag "$meta.id"
    label 'process_medium'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine != 'singularity' && !task.ext.singularity_pull_docker_container ?
    'biocontainers/bioconductor-updhmm:1.6.0--r44hdfd78af_0' :
    'https://depot.galaxyproject.org/singularity/bioconductor-updhmm:1.6.0--r44hdfd78af_0' }"

    input:
    tuple val(meta), path(upd_events_rds)

    output:
    tuple val(meta), path("*.upd_collapsed.txt")   , emit: upd_collapsed_txt
    tuple val(meta), path("*.upd_collapsed.rds")   , emit: upd_collapsed_rds
    tuple val(meta), path("*.R_sessionInfo.log")   , emit: session_info
    tuple val("${task.process}"), val("bioconductor-updhmm"), eval("Rscript -e \"library(UPDhmm); cat(as.character(packageVersion('UPDhmm')))\""), emit: versions_updhmm, topic: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    template 'collapseevents.R'
    
    stub:
    """
    touch ${meta.id}.upd_collapsed.txt
    touch ${meta.id}.upd_collapsed.rds
    touch ${meta.id}.R_sessionInfo.log
    """
}