process FASTQSCREEN_BUILDFROMINDEX {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastq-screen:0.15.3--pl5321hdfd78af_0':
        'biocontainers/fastq-screen:0.15.3--pl5321hdfd78af_0'}"

    input:
    val(indexes)    // Flattened list of meta, path pairs from INDEX_MODULE.out.index.collect()
    val(aligner)    // 'bwa', 'bowtie', 'bowtie2', or 'minimap2'

    output:
    path("FastQ_Screen_Genomes"), emit: database
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Default to bowtie2, since this is also the default of FastQ Screen itself
    aligner = aligner ?: "bowtie2"

    // Convert flat collected [meta,path,...] list into tuples [[meta,path],[meta,path],...]
    def genome_index_tuples = indexes.collate(2)

    // Map aligner to a representative index file extension (used to detect basename)
    def index_extensions = [
        'bwa': '.amb',              // BWA creates: .amb, .ann, .bwt, .pac, .sa
        'bowtie2': '.rev.1.bt2',    // Bowtie2 creates: .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, .rev.2.bt2
        'bowtie': '.1.ebwt*',        // Bowtie creates: .1.ebwt, .2.ebwt, .3.ebwt, .4.ebwt, .rev.1.ebwt, and .rev.2.ebwt or .ebwtl for large inputs
        'minimap2': '.mmi'          // Minimap2 creates: .mmi
    ]
    def extension_pattern = index_extensions[aligner]
    if (!extension_pattern) {
        error "Unsupported aligner: ${aligner}. Supported: bwa, bowtie, bowtie2, minimap2"
    }
    """
    GENOME_DIR="FastQ_Screen_Genomes"
    mkdir -p "\${GENOME_DIR}"

    # Add FastQ Screen config header
    echo "# FastQ Screen Configuration - ${aligner} indices" > "\${GENOME_DIR}/fastq_screen.conf"

    # Loop through genome, index directory pairs
    while read GENOME INDEX_DIR; do

        # copy index files to FastQ Screen database subdirectory
        OUTPUT_DIR="\$GENOME_DIR/\$GENOME"
        mkdir -p "\$OUTPUT_DIR"
        cp -r "\$INDEX_DIR/"* "\$OUTPUT_DIR/"

        #
        # Alternative approach that skips validation
        # Append to config using meta.id directly
        # echo -e "DATABASE\t\$GENOME\t\$OUTPUT_DIR" >> "\$GENOME_DIR/fastq_screen.conf"
        #

        # Find a representative index file to validate basename
        INDEX_FILE=\$(find "\$OUTPUT_DIR" -type f -name "*${extension_pattern}" | head -n1)
        if [ -z "\$INDEX_FILE" ]; then
            echo "ERROR: No ${aligner} index file (*${extension_pattern}) found in \$OUTPUT_DIR"
            exit 1
        fi

        # Extract basename
        if [ "${aligner}" = "bowtie" ]; then
            INDEX_BASE=\$(basename "\$INDEX_FILE" | sed 's/\\.1\\.ebwtl\\?\$//')
        else
            INDEX_BASE=\$(basename "\$INDEX_FILE" ${extension_pattern})
        fi

        # Validate that basename matches genome
        if [ "\$INDEX_BASE" != "\$GENOME" ]; then
            echo "ERROR: Index file basename (\$INDEX_BASE) does not match genome name (\$GENOME)"
            exit 1
        fi

        # Append genome and index path to FastQ Screen config
        echo -e "DATABASE\t\$GENOME\t\$OUTPUT_DIR/\$INDEX_BASE" >> "\$GENOME_DIR/fastq_screen.conf"

    done <<'EOF'
${genome_index_tuples.collect { meta, idx -> "${meta.id} ${idx}" }.join("\n")}
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqscreen: \$(echo \$(fastq_screen --version 2>&1) | sed 's/^.*FastQ Screen v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def genome_dir = "FastQ_Screen_Genomes"
    """
    mkdir ${genome_dir}
    touch ${genome_dir}/fastq_screen.conf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqscreen: \$(echo \$(fastq_screen --version 2>&1) | sed 's/^.*FastQ Screen v//; s/ .*\$//')
    END_VERSIONS
    """
}
