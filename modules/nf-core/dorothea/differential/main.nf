process DOROTHEA_DIFFERENTIAL_LIMMA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-decoupler:2.8.0--r43hdfd78af_0':
        'biocontainers/bioconductor-decoupler:2.8.0--r43hdfd78af_0'}"

    input:
    
    tuple val(meta), path(differential_result), path(pkn_network)

    output:
    
    tuple val(meta), path("*.dorothea.results.tsv")          , emit: results
    tuple val(meta), path("*.R_sessionInfo.log")          , emit: session_info
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    
    template 'dorothea_limma.R'
}
