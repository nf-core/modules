
process STADENIOLIB_CRAMFILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::staden_io_lib=1.14.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/staden_io_lib:1.14.14--h0d9da7e_3':
        'biocontainers/staden_io_lib:1.14.14--h0d9da7e_3' }"

    input:
    tuple val(meta), path(incram), path(cramIdx)
    val from
    val to

    output:
    tuple val(meta), path("*.cram"), emit: outcram
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cram_filter -n ${from}-${to} $args $incram ${prefix}_filtered.cram  

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stadeniolib: \$(cram_filter -h | head -n 1 |sed 's/^.*version //')
    END_VERSIONS
    """
}
