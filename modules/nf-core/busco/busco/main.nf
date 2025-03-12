process BUSCO_BUSCO {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c6/c607f319867d96a38c8502f751458aa78bbd18fe4c7c4fa6b9d8350e6ba11ebe/data'
        : 'community.wave.seqera.io/library/busco_sepp:f2dbc18a2f7a5b64'}"

    input:
    tuple val(meta), path(fasta, stageAs:'tmp_input/*')
    val mode                              // Required:    One of genome, proteins, or transcriptome
    val lineage                           // Required:    lineage for checking against, or "auto/auto_prok/auto_euk" for enabling auto-lineage
    path busco_lineages_path              // Recommended: BUSCO lineages file - downloads if not set
    path config_file                      // Optional:    BUSCO configuration file
    val clean_intermediates               // Optional:    Remove intermediate files

    output:
    tuple val(meta), path("*-busco.batch_summary.txt")                                        , emit: batch_summary
    tuple val(meta), path("short_summary.*.txt")                                              , emit: short_summaries_txt   , optional: true
    tuple val(meta), path("short_summary.*.json")                                             , emit: short_summaries_json  , optional: true
    tuple val(meta), path("*-busco/logs/busco.log")                                           , emit: log                   , optional: true
    tuple val(meta), path("*-busco/*/run_*/full_table.tsv")                                   , emit: full_table            , optional: true
    tuple val(meta), path("*-busco/*/run_*/missing_busco_list.tsv")                           , emit: missing_busco_list    , optional: true
    tuple val(meta), path("*-busco/*/run_*/single_copy_proteins.faa")                         , emit: single_copy_proteins  , optional: true
    tuple val(meta), path("*-busco/*/run_*/busco_sequences")                                  , emit: seq_dir               , optional: true
    tuple val(meta), path("*-busco/*/translated_proteins")                                    , emit: translated_dir        , optional: true
    tuple val(meta), path("*-busco")                                                          , emit: busco_dir
    tuple val(meta), path("busco_downloads/lineages/*")                                       , emit: downloaded_lineages   , optional: true
    tuple val(meta), path("*-busco/*/run_*/busco_sequences/single_copy_busco_sequences/*.faa"), emit: single_copy_faa       , optional: true
    tuple val(meta), path("*-busco/*/run_*/busco_sequences/single_copy_busco_sequences/*.fna"), emit: single_copy_fna       , optional: true

    path "versions.yml"                                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (mode !in ['genome', 'proteins', 'transcriptome']) {
        error("Mode must be one of 'genome', 'proteins', or 'transcriptome'.")
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}-${lineage}"
    def busco_config = config_file ? "--config ${config_file}" : ''
    def busco_lineage = lineage in ['auto', 'auto_prok', 'auto_euk']
        ? lineage.replaceFirst('auto', '--auto-lineage').replaceAll('_', '-')
        : "--lineage_dataset ${lineage}"
    def busco_lineage_dir = busco_lineages_path ? "--download_path ${busco_lineages_path}" : ''
    def intermediate_files = [
        './*-busco/*/auto_lineage',
        './*-busco/*/**/{miniprot,hmmer,.bbtools}_output',
        './*-busco/*/prodigal_output/predicted_genes/tmp/',
    ]
    def clean_cmd = clean_intermediates ? "rm -fr ${intermediate_files.join(' ')}" : ''
    """
    # Fix Augustus for Apptainer
    ENV_AUGUSTUS=/opt/conda/etc/conda/activate.d/augustus.sh
    set +u
    if [ -z "\${AUGUSTUS_CONFIG_PATH}" ] && [ -f "\${ENV_AUGUSTUS}" ]; then
        source "\${ENV_AUGUSTUS}"
    fi
    set -u

    # If the augustus config directory is not writable, then copy to writeable area
    if [ ! -w "\${AUGUSTUS_CONFIG_PATH}" ]; then
        # Create writable tmp directory for augustus
        AUG_CONF_DIR=\$( mktemp -d -p \$PWD )
        cp -r \$AUGUSTUS_CONFIG_PATH/* \$AUG_CONF_DIR
        export AUGUSTUS_CONFIG_PATH=\$AUG_CONF_DIR
        echo "New AUGUSTUS_CONFIG_PATH=\${AUGUSTUS_CONFIG_PATH}"
    fi

    # Ensure the input is uncompressed
    INPUT_SEQS=input_seqs
    mkdir "\$INPUT_SEQS"
    cd "\$INPUT_SEQS"
    for FASTA in ../tmp_input/*; do
        if [ "\${FASTA##*.}" == 'gz' ]; then
            gzip -cdf "\$FASTA" > \$( basename "\$FASTA" .gz )
        else
            ln -s "\$FASTA" .
        fi
    done
    cd ..

    busco \\
        --cpu ${task.cpus} \\
        --in "\$INPUT_SEQS" \\
        --out ${prefix}-busco \\
        --mode ${mode} \\
        ${busco_lineage} \\
        ${busco_lineage_dir} \\
        ${busco_config} \\
        ${args}

    # clean up
    rm -rf "\$INPUT_SEQS"
    ${clean_cmd}
    # find and remove broken symlinks from the cleanup
    find . -xtype l -delete

    # Move files to avoid staging/publishing issues
    mv ${prefix}-busco/batch_summary.txt ${prefix}-busco.batch_summary.txt
    mv ${prefix}-busco/*/short_summary.*.{json,txt} . || echo "Short summaries were not available: No genes were found."
    mv ${prefix}-busco/logs/busco.log ${prefix}-busco.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}-${lineage}"
    def fasta_name = files(fasta).first().name - '.gz'
    """
    touch ${prefix}-busco.batch_summary.txt
    mkdir -p ${prefix}-busco/${fasta_name}/run_${lineage}/busco_sequences

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
