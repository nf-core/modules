process VUEGEN {
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cd/cd83e0a51f0f2b0d4796e4c30cced46a1204f463a8707dad70ccd681e4ab588e/data'
:         'community.wave.seqera.io/library/python_vuegen_xz:23a454c211b5c866' }"

    input:
    val input_type
    path input_path
    val report_type

    output:
    path "*report", emit: output_folder
    tuple val("${task.process}"), val('vuegen'), eval('python -c "import vuegen; print(vuegen.__version__)"'), emit: versions_vuegen, topic: versions

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
            echo "Running VueGen with config file: ${input_path}"
            vuegen --config ${input_path} --report_type ${report_type} \$QUARTO_CHECK_FLAG ${args}
        elif [ "${input_type}" == "directory" ]; then
            echo "Running VueGen with directory: ${input_path}"
            vuegen --directory ${input_path} --report_type ${report_type} \$QUARTO_CHECK_FLAG ${args}
        fi
        """

    stub:
    """
        echo "STUB MODE: Creating a generic report directory"
        mkdir -p report
        touch report/report.txt
        """
}
