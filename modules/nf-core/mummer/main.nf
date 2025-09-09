process MUMMER {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mummer:3.23--pl5262h1b792b2_12' :
        'biocontainers/mummer:3.23--pl5262h1b792b2_12' }"

    input:
    tuple val(meta), path(ref), path(query)

    output:
    tuple val(meta), path("*.coords"), emit: coords
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed_ref = ref.getName().endsWith(".gz") ? true : false
    def fasta_name_ref = ref.getName().replace(".gz", "")

    def is_compressed_query = query.getName().endsWith(".gz") ? true : false
    def fasta_name_query = query.getName().replace(".gz", "")
    def VERSION = '3.23' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    if [ "$is_compressed_ref" == "true" ]; then
        gzip -c -d $ref > $fasta_name_ref
    fi
    if [ "$is_compressed_query" == "true" ]; then
        gzip -c -d $query > $fasta_name_query
    fi
    mummer \\
        $args \\
        $fasta_name_ref \\
        $fasta_name_query \\
        > ${prefix}.coords

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mummer: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '3.23' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.coords

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mummer: $VERSION
    END_VERSIONS
    """
}
