process CUSTOM_DUMPRUNPARAMS {
    label 'process_low'

    output:
    path("params_mqc.tsv"), emit: mqc_tsv

    script:
    run_params=params
    """
    echo "${run_params}" > run_params.txt
    table_header=\$( sed 's/,/\t/g' <( echo "Parameter,col1" ) )
    ## TODO: maybe make replacement of `:` only when surrounded by a-z1-9 ?
    ## EX: test_data:[sarscov2:[genome: OR https://raw <- need to test! Doesn't always work!
    run_params=\$( sed 's#^\\[##g;s#\\]\$##g' run_params.txt | sed 's#, #\\n#g;s#:#\t#g' | sed 's#^#--#g' )

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
