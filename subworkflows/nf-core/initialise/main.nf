// This is a reimplementation of the basic initialise methods in the nf-core template
// The role of this subworkflow is to check the input parameters and print help messages
// Use this to start your nf-core pipeline

include { paramsHelp; paramsSummaryLog; paramsSummaryMap; validateParameters } from 'plugin/nf-validation'

workflow INITIALISE {

    take:
        version         // bool
        help            // bool
        validate_params // bool
        logo            // bool

    main:

        // Get version string for use with --version or printing to log when creating
        String version_string = ""

        if (workflow.manifest.version) {
            def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
            version_string += "${prefix_v}${workflow.manifest.version}"
        }

        if (workflow.commitId) {
            def git_shortsha = workflow.commitId.substring(0, 7)
            version_string += "-g${git_shortsha}"
        }

        // Print workflow version and exit on --version
        if (version) {
            log.info "${workflow.manifest.name} ${version_string}"
            System.exit(0)
        }

        if ( logo ) {
            log.info logo(workflow, params.monochrome_logs)
        }
        log.info "${workflow.manifest.name} ${version_string}"

        // Print citation for nf-core
        def citation = "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  ${workflow.manifest.homePage}/blob/master/CITATIONS.md"
        log.info citation

        // Print help message if needed
        if (params.help) {
            def String command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv -profile docker"
            log.info paramsHelp(command)
            System.exit(0)
        }

        if ( params.validate_params != false ){
            validateParameters()
        }

        if (workflow.profile == 'standard' && workflow.configFiles.size() <= 1) {
            log.warn "[$workflow.manifest.name] You are attempting to run the pipeline without any custom configuration!\n\n" +
                    "This will be dependent on your local compute environment but can be achieved via one or more of the following:\n" +
                    "   (1) Using an existing pipeline profile e.g. `-profile docker` or `-profile singularity`\n" +
                    "   (2) Using an existing nf-core/configs for your Institution e.g. `-profile crick` or `-profile uppmax`\n" +
                    "   (3) Using your own local custom config e.g. `-c /path/to/your/custom.config`\n\n" +
                    "Please refer to the quick start section and usage docs for the pipeline.\n "
        }

        log.info paramsSummaryLog(workflow)

        summary_params = paramsSummaryMap(workflow)

    emit:
        summary_params

}

//
// ANSII Colours used for terminal logging
//
def logColours(Boolean monochrome_logs) {
    Map colorcodes = [:]

    // Reset / Meta
    colorcodes['reset']      = monochrome_logs ? '' : "\033[0m"
    colorcodes['bold']       = monochrome_logs ? '' : "\033[1m"
    colorcodes['dim']        = monochrome_logs ? '' : "\033[2m"
    colorcodes['underlined'] = monochrome_logs ? '' : "\033[4m"
    colorcodes['blink']      = monochrome_logs ? '' : "\033[5m"
    colorcodes['reverse']    = monochrome_logs ? '' : "\033[7m"
    colorcodes['hidden']     = monochrome_logs ? '' : "\033[8m"

    // Regular Colors
    colorcodes['black']      = monochrome_logs ? '' : "\033[0;30m"
    colorcodes['red']        = monochrome_logs ? '' : "\033[0;31m"
    colorcodes['green']      = monochrome_logs ? '' : "\033[0;32m"
    colorcodes['yellow']     = monochrome_logs ? '' : "\033[0;33m"
    colorcodes['blue']       = monochrome_logs ? '' : "\033[0;34m"
    colorcodes['purple']     = monochrome_logs ? '' : "\033[0;35m"
    colorcodes['cyan']       = monochrome_logs ? '' : "\033[0;36m"
    colorcodes['white']      = monochrome_logs ? '' : "\033[0;37m"

    // Bold
    colorcodes['bblack']     = monochrome_logs ? '' : "\033[1;30m"
    colorcodes['bred']       = monochrome_logs ? '' : "\033[1;31m"
    colorcodes['bgreen']     = monochrome_logs ? '' : "\033[1;32m"
    colorcodes['byellow']    = monochrome_logs ? '' : "\033[1;33m"
    colorcodes['bblue']      = monochrome_logs ? '' : "\033[1;34m"
    colorcodes['bpurple']    = monochrome_logs ? '' : "\033[1;35m"
    colorcodes['bcyan']      = monochrome_logs ? '' : "\033[1;36m"
    colorcodes['bwhite']     = monochrome_logs ? '' : "\033[1;37m"

    // Underline
    colorcodes['ublack']     = monochrome_logs ? '' : "\033[4;30m"
    colorcodes['ured']       = monochrome_logs ? '' : "\033[4;31m"
    colorcodes['ugreen']     = monochrome_logs ? '' : "\033[4;32m"
    colorcodes['uyellow']    = monochrome_logs ? '' : "\033[4;33m"
    colorcodes['ublue']      = monochrome_logs ? '' : "\033[4;34m"
    colorcodes['upurple']    = monochrome_logs ? '' : "\033[4;35m"
    colorcodes['ucyan']      = monochrome_logs ? '' : "\033[4;36m"
    colorcodes['uwhite']     = monochrome_logs ? '' : "\033[4;37m"

    // High Intensity
    colorcodes['iblack']     = monochrome_logs ? '' : "\033[0;90m"
    colorcodes['ired']       = monochrome_logs ? '' : "\033[0;91m"
    colorcodes['igreen']     = monochrome_logs ? '' : "\033[0;92m"
    colorcodes['iyellow']    = monochrome_logs ? '' : "\033[0;93m"
    colorcodes['iblue']      = monochrome_logs ? '' : "\033[0;94m"
    colorcodes['ipurple']    = monochrome_logs ? '' : "\033[0;95m"
    colorcodes['icyan']      = monochrome_logs ? '' : "\033[0;96m"
    colorcodes['iwhite']     = monochrome_logs ? '' : "\033[0;97m"

    // Bold High Intensity
    colorcodes['biblack']    = monochrome_logs ? '' : "\033[1;90m"
    colorcodes['bired']      = monochrome_logs ? '' : "\033[1;91m"
    colorcodes['bigreen']    = monochrome_logs ? '' : "\033[1;92m"
    colorcodes['biyellow']   = monochrome_logs ? '' : "\033[1;93m"
    colorcodes['biblue']     = monochrome_logs ? '' : "\033[1;94m"
    colorcodes['bipurple']   = monochrome_logs ? '' : "\033[1;95m"
    colorcodes['bicyan']     = monochrome_logs ? '' : "\033[1;96m"
    colorcodes['biwhite']    = monochrome_logs ? '' : "\033[1;97m"

    return colorcodes
}

//
// Does what is says on the tin
//
def dashedLine(monochrome_logs) {
    Map colors = logColours(monochrome_logs)
    return "-${colors.dim}----------------------------------------------------${colors.reset}-"
}

//
// nf-core logo
//
def logo(workflow, monochrome_logs) {
    Map colors = logColours(monochrome_logs)
    String.format(
        """\n
        ${dashedLine(monochrome_logs)}
                                                ${colors.green},--.${colors.black}/${colors.green},-.${colors.reset}
        ${colors.blue}        ___     __   __   __   ___     ${colors.green}/,-._.--~\'${colors.reset}
        ${colors.blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${colors.yellow}}  {${colors.reset}
        ${colors.blue}  | \\| |       \\__, \\__/ |  \\ |___     ${colors.green}\\`-._,-`-,${colors.reset}
                                                ${colors.green}`._,._,\'${colors.reset}
        ${dashedLine(monochrome_logs)}
        """.stripIndent()
    )
}
