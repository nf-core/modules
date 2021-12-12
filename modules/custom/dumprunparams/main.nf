process CUSTOM_DUMPRUNPARAMS {
    label 'process_low'

    input:
    val(exclude)

    output:
    path("params_mqc.tsv"), emit: mqc_tsv

    script:
    run_params = params
    to_exclude = exclude instanceof List && exclude.size() > 0 ? exclude.join('|') : ''
    exclude_cmd = exclude ? "| sed -E '/${to_exclude}/d'" : ''
    println exclude
    println exclude_cmd
    """
    echo "${run_params}" > run_params.txt
    table_header=\$( sed 's/,/\t/g' <( echo "Parameter,col1" ) )

    ## Remove outer brackets, remove nested arrays, split each param new line, delete lines containing params to ignore, convert params to cli flags
    run_params=\$( sed 's#^\\[##g;s#\\]\$##g' run_params.txt | sed 's/\\[.*]/\\.\\.\\.truncated_list\\.\\.\\./g' | sed 's#, #\\n#g;s#:#\t#g' ${exclude_cmd} | sed 's#^#--#g' )

    cat <<-MQC_HEADER >> params_mqc.tsv
    # plot_type: 'table'
    # section_name: 'Pipeline run parameters'
    # description: 'Resolved parameters for pipeline run: ${workflow.runName}. Note: All possible parameters are listed, but not necessarily used - some will only be utilised if a given module has been explicitly activated.'
    # pconfig:
    #     namespace: 'Cust Data'
    # headers:
    #     col1:
    #         title: 'Value'
    #         description: 'Final resolved values for all pipeline parameters'
    \$table_header
    \$run_params
    MQC_HEADER
    """
}
