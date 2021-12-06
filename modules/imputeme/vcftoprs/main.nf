process IMPUTEME_VCFTOPRS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "YOUR-TOOL-HERE" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/imputeme/vv1.0.7_cv1/imputeme_vv1.0.7_cv1.img' :
        'biocontainers/imputeme:vv1.0.7_cv1' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("id*/*.json"), emit: json
    path "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    #Set configuration - either from args or from defaults
    source("/imputeme/code/impute-me/functions.R")
    if(file.exists('$args')){
        set_conf("set_from_file",'$args')
    }else{
        set_conf("set_from_file", "/imputeme/code/impute-me/template/nextflow_default_configuration.R")
    }

    #main run
    return_message <- prepare_individual_genome('$vcf',overrule_vcf_checks=T)
    uniqueID <- sub(' </b>.+\$','',sub('^.+this run is <b> ','',return_message))
    convert_vcfs_to_simple_format(uniqueID=uniqueID)
    crawl_for_snps_to_analyze(uniqueIDs=uniqueID)
    run_export_script(uniqueIDs=uniqueID)

    #version export.
    f <- file("versions.yml","w")
    writeLines("${task.process}", f)
    writeLines(paste0(" imputeme: ", sub("^v","",get_conf("version"))),f)
    close(f)

    """

}
