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
        // container "quay.io/lassefolkersen/imputeme:v1.0.6"
        // TODO not fixed yet
    } else {
        // TODO - change to the biocontainer-based location
        container "quay.io/lassefolkersen/imputeme:latest"
    }

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.json"), emit: json
    path "*.version.txt"          , emit: version


    // The impute.me docker runs as user ubuntu, in /home/ubuntu. This works on direct AWS-implementations and on
    // handheld docker. However, when running in the nf-core setup (and also singularity fwiw) the home folder
    // is overwritten. The obvious fix, of using root or another user is a bit problematic because that strongly
    // affects the shiny-server setup which needs (read and write) access and complains about running as root.
    // For now, the most minimal solution seemed to be export of all required folders as volumes or tmpfs. In next
    // version of impute.me it's possible that the /home/ubuntu prefix will be made configurable, so many of these
    // exports can be avoided.
//    containerOptions "\
//        --mount 'type=tmpfs,source=,target=/home/ubuntu/logs' \
//        --mount 'type=volume,source=,target=/home/ubuntu/misc_files' \
//        --mount 'type=volume,source=,target=/home/ubuntu/configuration' \
//        --mount 'type=tmpfs,source=,target=/home/ubuntu/data' \
//        --mount 'type=volume,source=,target=/home/ubuntu/programs' \
//        --mount 'type=volume,source=,target=/home/ubuntu/prs_dir' \
//        --mount 'type=tmpfs,source=,target=/home/ubuntu/imputations' \
//        --mount 'type=tmpfs,source=,target=/home/ubuntu/vcfs' \
//        --mount 'type=volume,source=,target=/home/ubuntu/srv'"


    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "$options.args" variable
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    """
    #!/usr/bin/env Rscript

    #hacking test file to conform to expected test-results (shouldn't effect production running, and needed for short nf-core vcf files)
    #can be deleted if differently formed test-data is obtained
    to_insert<-"##fileformat=VCFv4.2"
    file.copy("$vcf","original.vcf.gz")
    system(paste0("zcat original.vcf.gz | sed '1i ",to_insert,"' | gzip -c > $vcf"))

    #set more verbose - this block re-writes the default configuration file to be more verbose
    #it's not really needed, other than for debugging, so this block can also be removed.
    system("echo 'verbose <- 10' > configuration.R")
    system("echo 'running_as_docker <- FALSE' >> configuration.R")
    system("echo 'block_double_uploads_by_md5sum <- FALSE' >> configuration.R")

    #setup minimal environment for vcf-processing
    #dir.create("~/logs/submission")
    source("/imputeme/code/impute-me/functions.R")

    #main run
    prepare_individual_genome('$vcf',overrule_vcf_checks=T,predefined_uniqueID="id_111111111")
    convert_vcfs_to_simple_format(uniqueID="id_111111111")
    crawl_for_snps_to_analyze(uniqueIDs="id_111111111")
    run_export_script(uniqueIDs="id_111111111")
    file.copy("data/id_111111111/id_111111111_data.json","output.json")

    #version export. Next impute-me software version has this as an internal function
    version_file_path="${software}.version.txt"
    f2<-file(version_file_path,"w")
    writeLines("v1.0.6",f2)
    close(f2)
    """
}
