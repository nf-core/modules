//
// Subworkflow that uses the nf-schema plugin to validate parameters and render the parameter summary
//

include { paramsSummaryLog   } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'

workflow UTILS_NFSCHEMA_PLUGIN {

    take:
    validate_params // boolean: validate the parameters

    main:

    //
    // Print parameter summary to stdout. This will display the parameters
    // that differ from the default given in the JSON schema
    //
    log.info paramsSummaryLog(workflow)

    //
    // Validate the parameters using nextflow_schema.json or the schema
    // given via the validation.parametersSchema configuration option
    //
    if(validate_params) {
        validateParameters()
    }

    emit:
    dummy_emit = true
}

