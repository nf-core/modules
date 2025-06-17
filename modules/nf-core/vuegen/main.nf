process VUEGEN {
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "dtu_biosustain_dsp/vuegen:v0.3.2-nextflow"

    input:
        val input_type
        path input_path
        val report_type

    output:
        path "*report", emit: output_folder
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        """
        # Validate quarto_check flag if using a conda environment
        if [[ "${task.conda}" != "null" ]]; then
            QUARTO_CHECK_FLAG="--quarto_checks"
        else
            QUARTO_CHECK_FLAG=""
        fi

        # Execute VueGen based on the input type
        if [ "${input_type}" == "config" ]; then
            echo "Running VueGen with config file: $input_path"
            vuegen --config $input_path --report_type $report_type \$QUARTO_CHECK_FLAG $args
        elif [ "${input_type}" == "directory" ]; then
            echo "Running VueGen with directory: $input_path"
            vuegen --directory $input_path --report_type $report_type \$QUARTO_CHECK_FLAG $args
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vuegen: \$( vuegen --version | sed -e "s/vuegen //g" )
        END_VERSIONS
        """

    stub:
        """
        echo "STUB MODE: Creating a generic report directory"
        mkdir -p report
        touch report/report.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vuegen: \$( vuegen --version | sed -e "s/vuegen //g" )
        END_VERSIONS
        """
}
