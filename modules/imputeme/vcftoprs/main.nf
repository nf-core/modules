// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/software
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

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
        container "quay.io/lassefolkersen/imputeme:latest"
    }

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml"            , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "$options.args" variable
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    """
    #!/usr/bin/env Rscript

    #Set temp configurations (these are for development (instant-git-pull), and will disappear)
    #set_conf("submission_logs_path","./")
    #set_conf("verbose",10)
    getwd()->a
    setwd("/imputeme/code/impute-me/")
    system("git pull")
    setwd(a)
    source("/imputeme/code/impute-me/functions.R")


    #Set configurations
    set_conf("defaults")
    set_conf("data_path","$TMPDIR/")
    set_conf("vcfs_path","$TMPDIR/")
    set_conf("autostart_supercronic",FALSE)
    set_conf("minimum_required_variant_in_vcf_count",1000)
    set_conf("modules_to_compute","ethnicity,AllDiseases") #remember to add AllDiseases and prs

    #main run
    d<-prepare_individual_genome('$vcf',overrule_vcf_checks=T)
    uniqueID<-sub(' </b>.+\$','',sub('^.+this run is <b> ','',d))
    convert_vcfs_to_simple_format(uniqueID=uniqueID)
    crawl_for_snps_to_analyze(uniqueIDs=uniqueID)
    run_export_script(uniqueIDs=uniqueID)
    setwd(a) #odd why this is needed, but seems like it is
    file.copy(paste0("$TMPDIR/",uniqueID,"/",uniqueID,"_data.json"),"output.json")

    #version export.
    version_file_path="versions.yml"
    f<-file(version_file_path,"w")
    writeLines("IMPUTEME_VCFTOPRS")
    writeLines(paste0(" imputeme: ", get_conf("version")),f)
    close(f)
    """


}
