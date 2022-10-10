process BBMAP_FILTERBYNAME {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bbmap=39.00" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.00--h5c4e2a8_0':
        'quay.io/biocontainers/bbmap:39.00--h5c4e2a8_0' }"

    input:
    tuple val(meta), path(reads)
    val output_extension
    val names

    output:
    tuple val(meta), path("*${output_extension}"), emit: reads
    tuple val(meta), path('*.log')               , emit: log
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def input = meta.single_end ? "in=${reads}" : "in=${reads[0]} in2=${reads[1]}"
    def filtered  = meta.single_end ? "out=${prefix}${output_extension}" : "out1=${prefix}_1${output_extension} out2=${prefix}_2${output_extension}"

    def avail_mem = 3
    if (!task.memory) {
        log.info '[filterbyname] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    filterbyname.sh \\
        names=$names \\
        -Xmx${avail_mem}g \\
        $input \\
        $filtered \\
        $args \\
        &> ${prefix}.filterbyname.log


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh)
    END_VERSIONS
    """
}
