process BBMAP_FILTERBYNAME {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.10--h92535d8_0':
        'biocontainers/bbmap:39.10--h92535d8_0' }"

    input:
    tuple val(meta), path(reads)
    val(names_to_filter)
    val(output_format)
    val(interleaved_output)

    output:
    tuple val(meta), path("*.${output_format}"), emit: reads
    tuple val(meta), path('*.log')             , emit: log
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input  = meta.single_end ? "in=${reads}" : "in=${reads[0]} in2=${reads[1]}"
    def output = (meta.single_end || interleaved_output) ?
        "out=${prefix}.${output_format}" :
        "out1=${prefix}_1.${output_format} out2=${prefix}_2.${output_format}"
    def names_command = names_to_filter ? "names=${names_to_filter}": ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[filterbyname] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    filterbyname.sh \\
        -Xmx${avail_mem}g \\
        $input \\
        $output \\
        $names_command \\
        $args \\
        | tee ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def filtered = (meta.single_end || interleaved_output) ?
        "echo '' | gzip > ${prefix}.${output_format}" :
        "echo '' | gzip >${prefix}_1.${output_format} ; echo '' | gzip >${prefix}_2.${output_format}"

    """
    $filtered
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """

}
