process RTGTOOLS_FORMAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rtg-tools:3.12.1--hdfd78af_0':
        'biocontainers/rtg-tools:3.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(input1), path(input2), path(sam_rg)

    output:
    tuple val(meta), path("*.sdf"), emit: sdf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Additional argument parsing
    // Currently a string, passed straight into the rtg cli
    // Want to seperate some args out into discrete args

    // Check if boolean is explicitly true. If not, set to false.
    // Don't use assert here
    def discard_empty_reads = task.ext.discard_empty_reads == true ? '--discard-empty-reads' : ''
    def no_gzip = task.ext.no_gzip == true ? '--no-gzip' : ''
    def reverse_complement = task.ext.reverse_complement == true ? '--reverse-complement' : ''

    // Should be a float. Default not specified
    def subsample = task.ext.subsample
    if (subsample) { assert subsample instanceof Float }

    // Should be an integer, default is 0.
    def seed = task.ext.seed
    if (seed) { assert seed instanceof Integer }

    // Should be an integer, default is 0.
    def trim_start_bases = task.ext.trim_start_bases
    if (trim_start_bases) { assert trim_start_bases instanceof Integer }

    def trim_end_bases = task.ext.trim_end_bases
    if (trim_end_bases) { assert trim_end_bases instanceof Integer }

    // Should be an integer, default is 0.
    def start_quality_threshold = task.ext.start_quality_threshold
    if (start_quality_threshold) { assert start_quality_threshold instanceof Integer }
    def end_quality_threshold = task.ext.end_quality_threshold
    if (end_quality_threshold) { assert end_quality_threshold instanceof Integer }

    // Should be a valid quality format
    // Allowed values are [sanger, solexa, illumina] (Default is sanger)
    def quality_format = task.ext.quality_format in [null, 'sanger', 'solexa', 'illumina'] ? task.ext.quality_format : error("Invalid value for task.ext.quality_format. Must be one of: sanger, solexa, illumina. If not provided, default is sanger.")

    // Should be an integer, default is 0.
    def min_read_length = task.ext.min_read_length
    if (min_read_length) { assert min_read_length instanceof Integer }

    // Should be an integer - could be linked in some way to the docker args/system resources
    // Default is available cores
    def threads = task.ext.threads
    if (threads) { assert threads instanceof Integer }

    // Keep existing args
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def single = meta.containsKey("single_end") ? meta.single_end : true

    def input = single ? "${input1}" : "--left ${input1} --right ${input2}"
    def rg = sam_rg ? "--sam-rg ${sam_rg}" : ""

    def avail_mem = "3G"
    if (!task.memory) {
        log.info '[RTG format] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue() + "M"
    }

    """
    rtg RTG_MEM=${avail_mem} format \\
        ${args} \\
        ${discard_empty_reads} \\
        ${no_gzip} \\
        ${reverse_complement} \\
        --quality-format=${quality_format} \\
        --min-read-length=${min_read_length} \\
        --threads=${threads} \\
        ${rg} \\
        --output ${prefix}.sdf \\
        ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rtg-tools: \$(echo \$(rtg version | head -n 1 | awk '{print \$4}'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = "3G"
    if (!task.memory) {
        log.info '[RTG format] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue() + "M"
    }
    """
    touch ${prefix}.sdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rtg-tools: \$(echo \$(rtg version | head -n 1 | awk '{print \$4}'))
    END_VERSIONS
    """
}
