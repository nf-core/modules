// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BBMAP_FILTERBYNAME {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c5/c5f02126e0d84b2eb0bcfd4b59c53fb192830fdf6510158c1a3cc121a7d3fff1/data'
        : 'community.wave.seqera.io/library/bbmap_pigz:7ce06264033efd7e'}"

    input:
    tuple val(meta), path(reads, stageAs: 'input/*'), path(names_file, stageAs: 'input/*')
    val(names_to_filter)
    val(output_format)
    val(interleaved_output)

    output:
    tuple val(meta), path(reads_out)      , emit: reads
    tuple val(meta), path("${prefix}.log"), emit: log
    tuple val("${task.process}"), val('bbmap'), eval('bbversion.sh | grep -v "Duplicate cpuset"'), emit: versions_bbmap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input  = meta.single_end ? "in=\"${reads}\"" : "in=\"${reads[0]}\" in2=\"${reads[1]}\""
    reads_out = (meta.single_end || interleaved_output) ?
        ["${prefix}.${output_format}"] :
        ["${prefix}_1.${output_format}", "${prefix}_2.${output_format}"]
    def output = reads_out.size() == 1 ?
        "out=\"${reads_out[0]}\"" :
        "out1=\"${reads_out[0]}\" out2=\"${reads_out[1]}\""
    // filterbyname.sh passes its arguments through eval, so a value with spaces needs its own inner quotes
    def names = [names_to_filter, names_file].flatten().findAll { name -> name }.join(',')
    def names_command = names ? "names='\"${names}\"'" : ""

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
        2>| >(tee "${prefix}.log" >&2)
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    reads_out = (meta.single_end || interleaved_output) ?
        ["${prefix}.${output_format}"] :
        ["${prefix}_1.${output_format}", "${prefix}_2.${output_format}"]
    def filtered = reads_out.collect { out_file ->
        output_format.endsWith('.gz')  ? "echo '' | gzip > ${out_file}" :
        output_format.endsWith('.bz2') ? "echo '' | bzip2 > ${out_file}" :
                                         "touch ${out_file}"
    }.join('\n')

    """
    $filtered
    touch ${prefix}.log
    """

}
