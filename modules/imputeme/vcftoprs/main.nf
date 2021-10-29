// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'


params.options = [:]
options        = initOptions(params.options)

process IMPUTEME_VCFTOPRS {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "YOUR-TOOL-HERE" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {

        // TODO not fixed yet - biocontainers only carry the docker image
    } else {
        // TODO - change to the biocontainer-based location (after code-changes to imputeme are done and up there)
        // container "quay.io/lassefolkersen/imputeme:latest"
        container "lassefolkersen/impute-me:latest"
    }

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml"            , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    #!/usr/bin/env Rscript

    #Set configuration - either from options.args or from defaults
    source("/imputeme/code/impute-me/functions.R")
    if(file.exists('$options.args')){
        set_conf("set_from_file",'$options.args')
    }else{
        set_conf("set_from_file", "/imputeme/code/impute-me/template/nextflow_default_configuration.R")
    }

    #main run
    return_message <- prepare_individual_genome('$vcf',overrule_vcf_checks=T)
    uniqueID <- sub(' </b>.+\$','',sub('^.+this run is <b> ','',return_message))
    convert_vcfs_to_simple_format(uniqueID=uniqueID)
    crawl_for_snps_to_analyze(uniqueIDs=uniqueID)
    run_export_script(uniqueIDs=uniqueID)
    file.copy(paste0("./",uniqueID,"/",uniqueID,"_data.json"),"output.json")

    #version export. Have to hardcode process name and software name because
    #won't run inside an R-block
    version_file_path="versions.yml"
    f <- file(version_file_path,"w")
    writeLines("IMPUTEME_VCFTOPRS:", f)
    writeLines(paste0(" imputeme: ", sub("^v","",get_conf("version"))),f)
    close(f)

    """

}
