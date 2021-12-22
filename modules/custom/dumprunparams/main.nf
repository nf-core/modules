include { dump_params_yml; indent_code_block } from "./parametrize"

def VERSION="0.0.1"

process CUSTOM_DUMPRUNPARAMS {
    label 'process_low'

    input:
    val(exclude)

    output:
    //path "params_mqc.yml", emit: mqc_yml
    //path "versions.yml"  , emit: versions

    // NOTE: removing script meant input channels were not detected
    script:
    run_params = params
    to_exclude = run_params.subMap(exclude)
    run_params_cleaned = exclude ? params.minus(to_exclude) : params
    run_params_formatted = run_params_cleaned.each{ it.getValue() instanceof LinkedHashMap | it.getValue() instanceof List | it.getValue() instanceof nextflow.config.ConfigMap ? run_params[it.getKey()] = '<truncated_list>' : null }
    println(run_params_formatted)

    // From https://stackoverflow.com/a/51288381/11502856
    // TODO try creating an importing custom groovy functions/methods (from link above) from a separate functions file


    // mqc_yml = [
    //     id: 'run_parameters',
    //     section_name: 'Pipeline run parameters',
    //     description: 'Resolved para,meters for pipeline run: ${workflow.runName}. Note: All possible parameters are listed, but not necessarily used - some will only be utilised if a given module has been explicitly activated.',
    //     plot_type: 'table',
    //     data: run_params_flattened
    // ]

    // WRITE FILE: https://www.nextflow.io/docs/latest/script.html#basic-read-write

    """
    """
}
