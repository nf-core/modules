process DOROTHEA_DIFFERENTIAL_LIMMA {
    tag "$meta.id"
    label 'dorothea_limma'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    
    tuple val(meta), path(differential_result)

    output:
    
    tuple val(meta), path("*.dorothea.results.tsv")          , emit: results
    tuple val(meta), path("*.dorothea.TF_plot.png")  , emit: TF_plot
    tuple val(meta), path("*.R_sessionInfo.log")          , emit: session_info
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    
    template 'dorothea_limma.R'
}
