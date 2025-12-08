process PLASTID_MAKE_WIGGLE {
    tag "$meta.id"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plastid:0.6.1--py39had3e4b6_2':
        'biocontainers/plastid:0.6.1--py39had3e4b6_2' }"

    input:
    tuple val(meta), path(bam), path(bam_index), path(p_offsets)
    val output_format

    output:
    tuple val(meta), path("*.{wig,bedgraph}"), emit: tracks
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = output_format == "bedgraph" ? "bedgraph" : "wig"
    def args = task.ext.args ?: ""
    """
    make_wiggle \\
        --count_files "$bam" \\
        --offset "$p_offsets" \\
        --fiveprime_variable \\
        --output_format "$output_format" \\
        -o "$prefix" \\
        $args

    if [ "$output_format" = "bedgraph" ]; then
        for FILE in *.wig; do
            mv "\$FILE" "\${FILE%.wig}.bedgraph"
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plastid: 0.6.1
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = output_format == "bedgraph" ? "bedgraph" : "wig"
    """
    touch ${prefix}_fw.$extension
    touch ${prefix}_rc.$extension

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plastid: 0.6.1
    END_VERSIONS
    """
}
