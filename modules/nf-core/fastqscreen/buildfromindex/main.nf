process FASTQSCREEN_BUILDFROMINDEX {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastq-screen:0.15.3--pl5321hdfd78af_0':
        'biocontainers/fastq-screen:0.15.3--pl5321hdfd78af_0'}"

    input:
    val(indexes)  // Flattened list of meta maps and index directory paths from INDEX_MODULE.out.index.collect()
    val(aligner)     // 'bwa', 'bowtie', 'bowtie2', or 'minimap2'

    output:
    path("FastQ_Screen_Genomes"), emit: database
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome_dir = "FastQ_Screen_Genomes"
    def index_tuples = indexes.collate(2)    // create pairs of meta, dir tuples
    def genome_names = index_tuples.collect { meta, _idx -> meta.id }
    def index_paths = index_tuples.collect { _meta, idx -> idx }

    // Default to bowtie2, since this is also the default of FastQ Screen itself
    aligner = aligner ?: "bowtie2"

    // Map aligner to a representative index file extension (used to detect basename)
    def index_extensions = [
        'bwa': '.amb',              // BWA creates: .amb, .ann, .bwt, .pac, .sa
        'bowtie2': '.rev.1.bt2',    // Bowtie2 creates: .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, .rev.2.bt2
        'bowtie': '1.ebwt*',        // Bowtie creates: .1.ebwt, .2.ebwt, .3.ebwt, .4.ebwt, .rev.1.ebwt, and .rev.2.ebwt or .ebwtl for large inputs
        'minimap2': '.mmi'          // Minimap2 creates: .mmi
    ]
    def extension_pattern = index_extensions[aligner]
    if (!extension_pattern) {
        error "Unsupported aligner: ${aligner}. Supported: bwa, bowtie, bowtie2, minimap2"
    }
    """
    mkdir $genome_dir

    # Build config header
    echo "# FastQ Screen Configuration - ${aligner} indices" > ${genome_dir}/fastq_screen.conf

    # Create arrays of genome names and index paths
    GENOMES=(${genome_names.collect { "\"${it}\"" }.join(' ')})
    INDICES=(${index_paths.collect { "\"${it}\"" }.join(' ')})

    # Loop over genome and indices array using the array index to keep them in sync
    for i in "\${!GENOMES[@]}"; do

        # Create output directory named after the current genome
        GENOME="\${GENOMES[\$i]}"
        OUTPUT_DIR="${genome_dir}/\${GENOME}"

        # Copy all index files
        cp -r "\${INDICES[\$i]}" "\${OUTPUT_DIR}"

        # Use a representative file to detect the index basename (should match genome name)
        # All index files for a genome should share the same basename prefix
        INDEX_FILE=\$(find "\${OUTPUT_DIR}" -type f -name "*${extension_pattern}" | head -n1)

        if [ -z "\${INDEX_FILE}" ]; then
            echo "ERROR: No ${aligner} index file (*${extension_pattern}) found in \${OUTPUT_DIR}"
            exit 1
        fi

        # Extract basename by removing the extension
        # Note: basename approach does not work for bowtie where there
        # could be two different extensions based on the size of the genome.
        # Because of the potential presence of dots in genome names, we cannot use
        # \${INDEX_BASE%%.*} to remove everything after the first dot either.
        if [ ${aligner} == "bowtie" ]; then
            INDEX_BASE=\$(basename "\${INDEX_FILE}" | sed 's/\\.1\\.ebwtl\\?\$//')
        else
            INDEX_BASE=\$(basename "\${INDEX_FILE}" ${extension_pattern})
        fi

        # Check if basename matches genome name
        if [ \${INDEX_BASE} != \${GENOME} ]; then
            echo "ERROR: Filename of index files (\${INDEX_FILE} - \${INDEX_BASE}) does not match the expected genome name (\${GENOME})."
            exit 1
        fi

        # FastQ Screen will find all related index files using this basename in the config
        echo "DATABASE\t\${GENOME}\t\${OUTPUT_DIR}/\${INDEX_BASE}" >> ${genome_dir}/fastq_screen.conf
    done

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
