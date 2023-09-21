// This is a reimplementation of the basic initialise methods in the nf-core template
// The role of this subworkflow is to check the input parameters and print help messages
// Use this to start your nf-core pipeline

include { paramsHelp; paramsSummaryLog; validateParameters } from 'plugin/nf-validation'

workflow INITIALISE {

    // Print workflow version and exit on --version
    if (params.version) {
        String version_string = ""

        if (workflow.manifest.version) {
            def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
            version_string += "${prefix_v}${workflow.manifest.version}"
        }

        if (workflow.commitId) {
            def git_shortsha = workflow.commitId.substring(0, 7)
            version_string += "-g${git_shortsha}"
        }
        log.info "${workflow.manifest.name} ${version_string}"
        System.exit(0)
    }

    // Print help message if needed
    if (params.help) {
        def String command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv -profile docker"
        log.info paramsHelp(command)
        System.exit(0)
    }

    if ( params.validate_parameters != false ){
        validateParameters()
    }

    // Print citation for nf-core
    def citation = "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
        "* The nf-core framework\n" +
        "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
        "* Software dependencies\n" +
        "  ${workflow.manifest.homePage}/blob/master/CITATIONS.md"
    log.info citation

    if (workflow.profile == 'standard' && workflow.configFiles.size() <= 1) {
        log.warn "[$workflow.manifest.name] You are attempting to run the pipeline without any custom configuration!\n\n" +
                "This will be dependent on your local compute environment but can be achieved via one or more of the following:\n" +
                "   (1) Using an existing pipeline profile e.g. `-profile docker` or `-profile singularity`\n" +
                "   (2) Using an existing nf-core/configs for your Institution e.g. `-profile crick` or `-profile uppmax`\n" +
                "   (3) Using your own local custom config e.g. `-c /path/to/your/custom.config`\n\n" +
                "Please refer to the quick start section and usage docs for the pipeline.\n "
    }

    log.info paramsSummaryLog(workflow)
}
