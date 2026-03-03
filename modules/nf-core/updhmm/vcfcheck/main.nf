process UPDHMM_VCFCHECK {
    tag "$meta.id"
    label 'process_medium'
   
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine != 'singularity' && !task.ext.singularity_pull_docker_container ?
    'biocontainers/bioconductor-updhmm:1.6.0--r44hdfd78af_0' :
    'https://depot.galaxyproject.org/singularity/bioconductor-updhmm:1.6.0--r44hdfd78af_0' }"
   
    input:
    tuple val(meta), path(vcf), path(tbi)
 
    output:
    tuple val(meta), path("*.processed.rds")    , emit: processed_vcf
    tuple val(meta), path("*.R_sessionInfo.log"), emit: session_info
    tuple val("${task.process}"), val("bioconductor-updhmm"), eval("Rscript -e \"library(UPDhmm); cat(as.character(packageVersion('UPDhmm')))\""), emit: versions_updhmm, topic: versions
    tuple val("${task.process}"), val("bioconductor-variantannotation"), eval("Rscript -e \"library(VariantAnnotation); cat(as.character(packageVersion('VariantAnnotation')))\""), emit: versions_variantannotation, topic: versions
    tuple val("${task.process}"), val("r-base"), eval("Rscript -e \"cat(strsplit(version[['version.string']], ' ')[[1]][3])\""), emit: versions_rbase, topic: versions
   
    when:
    task.ext.when == null || task.ext.when
 
    script:
    template 'vcfcheck.R'
 
    stub:
    """
    touch ${meta.id}.processed.rds
    touch ${meta.id}.R_sessionInfo.log
    """
}
