process BUSCO {
    tag "$meta.id"
    label 'process_medium'
    
    conda (params.enable_conda ? "bioconda::busco=5.2.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.2.2--pyhdfd78af_0':
        'quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    val(mode)
    path(augustus_config)
    val(lineage)

    output:
    tuple val(meta), path("${meta.id}/run_*/full_table.tsv"), emit: tsv
    tuple val(meta), path("${meta.id}/run_*/short_summary.txt"), emit: txt
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (lineage) args += " --lineage_dataset $lineage"
    """
    # Ensure the input is uncompressed
    gzip -cdf $fasta > __UNCOMPRESSED_FASTA_FILE__
    
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
    
    # Copy the image's AUGUSTUS config directory if it was not provided to the module
    [ ! -e augustus_config ] && cp -a /usr/local/config augustus_config
    
    # Busco command
    AUGUSTUS_CONFIG_PATH=augustus_config busco $args --augustus --mode $mode --cpu $task.cpus --in __UNCOMPRESSED_FASTA_FILE__ --out $meta.id

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
