process BBMAP_FILTERBYNAME {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bbmap=39.01" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.01--h311275f_0':
        'quay.io/biocontainers/bbmap:39.01--h311275f_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.$output_extension"), emit: reads
    tuple val(meta), path('*.log')              , emit: log
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args   ?: [ args:"" ]
    args.args = args.args ? "$args.args" : ""

    def prefix             = task.ext.prefix ?: "${meta.id}"
    def interleaved_output = args.interleaved_output ? true : false
    def input  = meta.single_end ? "in=${reads}" :
        "in=${reads[0]} in2=${reads[1]}"

    def avail_mem = 3
    if (!task.memory) {
        log.info '[filterbyname] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    output_extension   = args.output_extension ? "$args.output_extension" :
        meta.single_end ? reads.name.tokenize('.')[1..-1].join('.') :
        reads[0].name.tokenize('.')[1..-1].join('.')

    filtered = (meta.single_end || interleaved_output) ?
        "out=${prefix}.$output_extension" :
        "out1=${prefix}_1.$output_extension out2=${prefix}_2.$output_extension"

    """
    filterbyname.sh \\
        -Xmx${avail_mem}g \\
        $input \\
        $filtered \\
        $args.args \\
        &> ${prefix}.filterbyname.log


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh)
    END_VERSIONS
    """
}
