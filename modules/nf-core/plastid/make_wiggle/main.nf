process PLASTID_MAKE_WIGGLE {
    tag "$meta.id"
    label "process_single"

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plastid:0.6.1--py39had3e4b6_2':
        'biocontainers/plastid:0.6.1--py39had3e4b6_2' }"

    input:
    tuple val(meta), path(bam), path(bam_index), path(p_offsets)
    val(mapping_rule)

    output:
    tuple val(meta), path("*.{wig,bedgraph}"), emit: tracks
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (mapping_rule == 'fiveprime_variable' && !p_offsets) {
        error "p_offsets file is required when using mapping_rule 'fiveprime_variable'"
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""
    def offset_arg = mapping_rule == 'fiveprime_variable' ? "--offset $p_offsets" : ""
    def extension = args.contains('--output_format bedgraph') ? "bedgraph" : "wig"
    def VERSION = "0.6.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    make_wiggle \\
        --count_files "$bam" \\
        $offset_arg \\
        --${mapping_rule} \\
        -o "$prefix" \\
        $args

    if [ "$extension" = "bedgraph" ]; then
        for FILE in *.wig; do
            mv "\$FILE" "\${FILE%.wig}.bedgraph"
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plastid: $VERSION
    END_VERSIONS
    """

    stub:
    if (mapping_rule == 'fiveprime_variable' && !p_offsets) {
        error "p_offsets file is required when using mapping_rule 'fiveprime_variable'"
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""
    def extension = args.contains('--output_format bedgraph') ? "bedgraph" : "wig"
    def VERSION = "0.6.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}_fw.${extension}
    touch ${prefix}_rc.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plastid: $VERSION
    END_VERSIONS
    """
}
