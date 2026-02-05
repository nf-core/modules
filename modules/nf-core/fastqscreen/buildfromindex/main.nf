process FASTQSCREEN_BUILDFROMINDEX {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastq-screen:0.15.3--pl5321hdfd78af_0':
        'biocontainers/fastq-screen:0.15.3--pl5321hdfd78af_0'}"

    input:
    val(genome_names)                       // Flattened list of [index_dir] from INDEX_MODULE.out.index.collect { meta, _index_dir -> meta.id }
    path(index_paths, stageAs: "index_?/*") // Flattened list of [index_dir] from INDEX_MODULE.out.index.collect { _meta, index_dir -> index_dir }
    val(aligner)     // 'bwa', 'bowtie', 'bowtie2', or 'minimap2'

    output:
    path("FastQ_Screen_Genomes"), emit: database
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Default to bowtie2, since this is also the default of FastQ Screen itself
    aligner = aligner ?: "bowtie2"

    def genome_dir = "FastQ_Screen_Genomes"

    // Validation: Check for duplicate genome names
    def unique_genomes = genome_names.unique()
    if (unique_genomes.size() != genome_names.size()) {
        def duplicates = genome_names.findAll { name ->
            genome_names.count(name) > 1
        }.unique()
        error "Duplicate genome names detected: ${duplicates.join(', ')}. Each genome must have a unique name."
    }

    // Validation: Check that number of genomes matches number of index paths
    if (genome_names.size() != index_paths.size()) {
        error "Mismatch: ${genome_names.size()} genome names provided but ${index_paths.size()} index paths. They must match."
    }

    // convert to space-separated input for bash; avoids literal list values to be passed (e.g. `"[GRCh38.chr21.fa," "PlasmoDB-68_Pfalciparum3D7_Genome]"`)
    genome_names = genome_names.collect { "\"$it\"" }.join(' ')

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
    mkdir ${genome_dir}

    # Build config header
    echo "# FastQ Screen Configuration - ${aligner} indices" > ${genome_dir}/fastq_screen.conf

    # Copy all index files to the same directory
    for idx in ${index_paths}; do
        cp "\${idx}/"* ${genome_dir}
    done

    # Create database lines in fastqscreen.conf based on supplied genome names
    for genome in ${genome_names}; do
        echo "DATABASE\t\${genome}\t${genome_dir}/\${genome}" >> ${genome_dir}/fastq_screen.conf
    done

    # Perform sanity checks
    # Note: these checks only check for the presence of correctly named files, but does
    # not guarantee that they are matched appropriately with their genomes
    for genome in ${genome_names}; do

        index_file=\$(find "${genome_dir}" -type f -name "\${genome}${extension_pattern}" | head -n1)

        # Check if the expected index files are present based on the chosen aligner
        # Use a representative file to detect the index basename (should match genome name)
        # All index files for a genome should share the same basename prefix
        if [ -z "\${index_file}" ]; then
            echo "ERROR: No ${aligner} index file (\${genome}${extension_pattern}) found in ${genome_dir}."
            exit 1
        fi

        # Check if the basename of the index files matches the genome basename
        # by removing the extension
        # Note: basename approach does not work for bowtie where there
        # could be two different extensions based on the size of the genome.
        # Because of the potential presence of dots in genome names, we cannot use
        # \${index_basename%%.*} to remove everything after the first dot either.
        if [ "${aligner}" == "bowtie" ]; then
            index_basename=\$(basename "\${index_file}" | sed 's/\\.1\\.ebwtl\\?\$//')
        else
            index_basename=\$(basename "\${index_file}" ${extension_pattern})
        fi
        if [ "\${index_basename}" != "\${genome}" ]; then
            echo "ERROR: Filename of index files (\${index_file} - \${index_basename}) does not match the expected genome name (\${genome})."
            exit 1
        fi
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
