process UPDHMM_CALCULATEEVENTS {
    tag "$meta.id"
    label 'process_high'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine != 'singularity' && !task.ext.singularity_pull_docker_container ?
    'biocontainers/bioconductor-updhmm:1.6.0--r44hdfd78af_0' :
    'https://depot.galaxyproject.org/singularity/bioconductor-updhmm:1.6.0--r44hdfd78af_0' }"

    input:
    tuple val(meta), path(processed_vcf_rds)

    output:
    tuple val(meta), path("*.upd_events.txt")   , emit: upd_events_txt
    tuple val(meta), path("*.upd_events.rds")   , emit: upd_events_rds
    tuple val(meta), path("*.R_sessionInfo.log"), emit: session_info
    
    tuple val("${task.process}"), val("bioconductor-updhmm"), eval("Rscript -e \"library(UPDhmm); cat(as.character(packageVersion('UPDhmm')))\""), emit: versions_updhmm, topic: versions
    tuple val("${task.process}"), val("bioconductor-biocparallel"), eval("Rscript -e \"library(BiocParallel); cat(as.character(packageVersion('BiocParallel')))\""), emit: versions_biocparallel, topic: versions
    tuple val("${task.process}"), val("r-base"), eval("Rscript -e \"cat(strsplit(version[['version.string']], ' ')[[1]][3])\""), emit: versions_rbase, topic: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    template 'calculateevents.R'
    
    stub:
    """
    touch ${meta.id}.upd_events.txt
    touch ${meta.id}.upd_events.rds
    touch ${meta.id}.R_sessionInfo.log
    """
}
