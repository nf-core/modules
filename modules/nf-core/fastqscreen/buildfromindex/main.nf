process FASTQSCREEN_BUILDFROMINDEX {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/19/19cd241b503facdddc020294a0996a5c8b2570f19ea57b3f404618ac5988a87f/data':
        'community.wave.seqera.io/library/bowtie2_bowtie_bwa_fastq-screen_minimap2:db243c473e0c4a39'}"

    input:
    path(index_paths, stageAs: "index_?/*")  // Flattened list of [index_dir] from INDEX_MODULE.out.index.collect{ _meta, dir -> dir }
    val(aligner)    // 'bwa', 'bowtie', 'bowtie2', or 'minimap2'

    output:
    path("FastQ_Screen_Genomes"), emit: database
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Default to bowtie2, since this is also the default of FastQ Screen itself
    aligner = aligner ?: "bowtie2"

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

    for INDEX_DIR in ${index_paths}; do

        # Validate that at least one index file exists given the specified aligner
        # All index files for a genome should share the same basename prefix
        # Note: since we are searching in the original input staged files (which are symlinks)
        # the -L flag is required
        INDEX_FILE=\$(find -L "\${INDEX_DIR}" -type f -name "*${extension_pattern}" | head -n1)
        if [ -z "\$INDEX_FILE" ]; then
            echo "ERROR: No ${aligner} index file matching ${extension_pattern} found in \${INDEX_DIR}."
            exit 1
        fi

        # Use the representative file to detect the genome name
        if [ "${aligner}" == "bowtie" ]; then
            GENOME=\$(basename "\${INDEX_FILE}" | sed 's/\\.1\\.ebwtl\\?\$//')
        else
            GENOME=\$(basename "\${INDEX_FILE}" ${extension_pattern})
        fi

        # Create output directory
        OUTPUT_DIR="\${GENOME_DIR}/\${GENOME}"
        mkdir -p "\${OUTPUT_DIR}"

        # Copy index files into the output directory
        cp -r "\${INDEX_DIR}/"* "\$OUTPUT_DIR/"

        # Append genome and index path to FastQ Screen config
        echo -e "DATABASE\t\${GENOME}\t\${OUTPUT_DIR}/\${GENOME}" >> "\${GENOME_DIR}/fastq_screen.conf"
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
