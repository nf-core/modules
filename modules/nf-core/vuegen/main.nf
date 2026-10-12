process VUEGEN {
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a0ff4f778cefa7ae78c684ba2b97dc4a1fac3d49ffc958c532bfe06e19e23807/data'
:         'community.wave.seqera.io/library/python_pytinytex_quarto_r-tinytex_pruned:9eebdec0448f6563' }"

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

        # Same Quarto cache fix as modules/nf-core/quartonotebook.
        # Fix Quarto for Apptainer (see https://community.seqera.io/t/confusion-over-why-a-tool-works-in-docker-but-fails-in-singularity-when-the-installation-doesnt-differ-i-e-using-wave-micromamba/1244)
        # HOME is also made writable so Quarto can install TinyTeX for PDF reports.
        # Same fix as in modules/nf-core/memote/report and modules/nf-core/bcftools/plotvcfstats.
        mkdir -p nxf_home
        export HOME=\$PWD/nxf_home
        export XDG_CACHE_HOME="./.xdg_cache_home"
        export XDG_DATA_HOME="./.xdg_data_home"
        ENV_QUARTO=/opt/conda/etc/conda/activate.d/quarto.sh
        set +u
        if [ -z "\${QUARTO_DENO}" ] && [ -f "\${ENV_QUARTO}" ]; then
            source "\${ENV_QUARTO}"
        fi
        set -u

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
