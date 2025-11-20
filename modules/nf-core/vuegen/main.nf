process VUEGEN {
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fa/fadd4c6459b24fc3964d47d72dbf809e425054e08f1aec9d56c8bec40b4b3a47/data'
    : 'community.wave.seqera.io/library/vuegen_python:236414fc5cfce774'}"

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
        # Set environment variables needed for Quarto rendering
        # (needed for apptainer/singularity)
        export XDG_CACHE_HOME="./.xdg_cache_home"
        export XDG_DATA_HOME="./.xdg_data_home"
        # Fix Quarto for apptainer: activate conda environment
        # https://github.com/mahesh-panchal/quarto-docker-singularity-problem
        ENV_QUARTO="\${ENV_QUARTO:-/opt/conda/etc/conda/activate.d/quarto.sh}"
        set +u
        if [ -z "\${QUARTO_DENO}" ] && [ -f "\${ENV_QUARTO}" ]; then
            source "\${ENV_QUARTO}"
        fi
        set -u

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
            vuegen: \$( python -c "import vuegen; print(vuegen.__version__)" )
        END_VERSIONS
        """

    stub:
        """
        echo "STUB MODE: Creating a generic report directory"
        mkdir -p report
        touch report/report.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vuegen: \$( python -c "import vuegen; print(vuegen.__version__)" )
        END_VERSIONS
        """
}
