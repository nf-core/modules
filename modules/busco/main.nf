process BUSCO {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::busco=5.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.3.2--pyhdfd78af_0':
        'quay.io/biocontainers/busco:5.3.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)  // Required:    meta map, and fasta sequence file
    each lineage                  // Required:    lineage to check against
    path busco_lineages_path      // Recommended: path to busco lineages - downloads if not set
    path config_file              // Optional:    busco configuration file

    output:
    tuple val(meta), path("*-busco"), emit: busco_dir
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}-${lineage}"
    def busco_config = config_file ? "--config $config_file" : ''
    def busco_lineage_dir = busco_lineages_path ? "--download_path ${busco_lineages_path}" : ''
    """
    # Nextflow changes the container --entrypoint to /bin/bash (container default entrypoint: /usr/local/env-execute)
    # Check for container variable initialisation script and source it.
    if [ -f "/usr/local/env-activate.sh" ]; then
        # . "/usr/local/env-activate.sh"  # Errors out because of various unbound variables
        export PATH='/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin'
        export CONDA_PREFIX='/usr/local'
        export CONDA_SHLVL='1'
        export CONDA_DEFAULT_ENV='/usr/local'
        export CONDA_PROMPT_MODIFIER=''
        . "/usr/local/etc/conda/activate.d/activate-r-base.sh"
        . "/usr/local/etc/conda/activate.d/augustus.sh"
        . "/usr/local/etc/conda/activate.d/openjdk_activate.sh"
    fi

    # If the augustus config directory is not writable, then copy to writeable area
    if [ ! -w "\${AUGUSTUS_CONFIG_PATH}" ]; then
        # Create writable tmp directory for augustus
        AUG_CONF_DIR=\$( mktemp -d -p \$PWD )
        cp -r \$AUGUSTUS_CONFIG_PATH/* \$AUG_CONF_DIR
        export AUGUSTUS_CONFIG_PATH=\$AUG_CONF_DIR
        echo "New AUGUSTUS_CONFIG_PATH=\${AUGUSTUS_CONFIG_PATH}"
    fi

    # Ensure the input is uncompressed
    gzip -cdf $fasta > ${prefix}_uncompressed.fasta

    busco \\
        --cpu $task.cpus \\
        --in ${prefix}_uncompressed.fasta \\
        --out ${prefix}-busco \\
        --lineage_dataset $lineage \\
        $busco_lineage_dir \\
        $busco_config \\
        $args

    # clean up
    rm ${prefix}_uncompressed.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
