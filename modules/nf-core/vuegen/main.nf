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
        # Validate Python version if using a conda environment
        if [[ "${task.conda}" != "null" ]]; then
            echo "Checking Python version in Conda environment..."
            PYTHON_VERSION=\$(python --version 2>&1 | awk '{print \$2}')
            echo "Python version: \$PYTHON_VERSION"
            REQUIRED_VERSION="3.11.0"

            # Check if Python version is lower than the required one
            if [ "\$(printf '%s\\n' \"\$PYTHON_VERSION\" \"\$REQUIRED_VERSION\" | sort -V | head -n1)" != "\$REQUIRED_VERSION" ]; then
                echo "ERROR: This module requires Python >= \$REQUIRED_VERSION. Current version: \$PYTHON_VERSION" >&2
                exit 1
            fi
        fi

        # Execute VueGen based on the input type
        if [ "${input_type}" == "config" ]; then
            echo "Running VueGen with config file: $input_path"
            vuegen --config $input_path --report_type $report_type $args
        elif [ "${input_type}" == "directory" ]; then
            echo "Running VueGen with directory: $input_path"
            vuegen --directory $input_path --report_type $report_type $args
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
