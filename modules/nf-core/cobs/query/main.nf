process COBS_QUERY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cobs:0.3.0--hdcf5f25_1' :
        'biocontainers/cobs:0.3.0--hdcf5f25_1'}"

    input:
    tuple val(query_meta), path(query)
    tuple val(index_meta), path(index)

    output:
    tuple val(meta), path("*.matches.gz"), emit: matches
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Note: xz is explicitly supported in this module as COBS indexes compress well under xz (see https://doi.org/10.1101/2023.04.15.536996)

    // variables setup
    meta = [id: query_meta.id, index_id: index_meta.id]
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${query_meta.id}-${index_meta.id}"
    def should_load_the_whole_index_into_RAM = args.contains("--load-complete")
    def index_is_gzip_compressed = index.toString().endsWith(".gz")
    def index_is_xz_compressed = index.toString().endsWith(".xz")
    def is_compact_and_compressed_index = index.toString().endsWith(".cobs_compact.gz") || index.toString().endsWith(".cobs_compact.xz")
    def decompress_tool = ""
    def get_index_size_command = ""
    if (index_is_gzip_compressed) {
        decompress_tool = "gzip"
        // note: unfortunately "--list" is not available in BusyBox gzip, so we need to decompress the archive entirely to get the index size
        // note: this is not ideal and should be optimised in next versions of this module, although I don't know how yet
        get_index_size_command = "index_size=\$(gzip -d -c $index | wc -c)"
    } else if (index_is_xz_compressed) {
        decompress_tool = "xz"
        get_index_size_command = "index_size=\$(xz --list --robot $index | grep file | awk '{print \$5}')"
    }
    def command = ""

    // input checking
    if (is_compact_and_compressed_index && should_load_the_whole_index_into_RAM) {
        error "COBS does not support loading a compact and compressed index directly RAM - please remove --load-complete from the args or use a classic index or decompress the index"
    }


    // run COBS
    def decompressed_index = ""
    def should_delete_decompressed_index = false
    if (index_is_gzip_compressed || index_is_xz_compressed) {
        if (should_load_the_whole_index_into_RAM) {
            // streams compressed index to RAM
            command += get_index_size_command
            command += """
                        cobs \\
                            query \\
                            $args \\
                            -T $task.cpus \\
                            -i <(${decompress_tool} -d -c "${index}") \\
                            -f <(zcat $query) \\
                            --index-sizes \$index_size \\
                       """
        } else {
            // decompresses compressed index to disk and mmap it
            decompressed_index = index.toString().replaceAll(".[gx]z\$", "")
            should_delete_decompressed_index = true
            command += """
                            ${decompress_tool} -d -c "${index}" > "${decompressed_index}"

                            cobs \\
                                query \\
                                $args \\
                                -T $task.cpus \\
                                -i $decompressed_index \\
                                -f <(zcat $query) \\
                       """
        }
    }else {
        // index is not compressed - will be loaded completely into RAM (if --load-complete is passed) or memory mapped otherwise
        command += """
                        cobs \\
                            query \\
                            $args \\
                            -T $task.cpus \\
                            -i $index \\
                            -f <(zcat $query) \\
                   """
    }
    command += " | gzip > ${prefix}.matches.gz"

    if (should_delete_decompressed_index) {
        command += """
                        rm -v "${decompressed_index}"
                   """
    }

    // note: for some reason the indentation is not stripped in some tests making the versions.yml file invalid
    // note: thus I am stripping it directly here
    command += """
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    cobs: \$(cobs version 2>&1 | awk '{print \$3}')
END_VERSIONS
    """
    command

    stub:
    meta = [id: query_meta.id, index_id: index_meta.id]
    def prefix = task.ext.prefix ?: "${query_meta.id}-${index_meta.id}"
    """
    # Note: just using touch yields java.io.EOFException when snapshotting the stub output - need to create a valid gzip file
    gzip </dev/null >${prefix}.matches.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobs: \$(cobs version 2>&1 | awk '{print \$3}')
    END_VERSIONS
    """
}
