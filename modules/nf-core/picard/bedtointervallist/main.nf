process PICARD_BEDTOINTERVALLIST {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::picard=2.27.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.27.4--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.27.4--hdfd78af_0' }"

    input:
    tuple val(meta) , path(bed)
    tuple val(meta2), path(dict)
    file  arguments_file

    output:
    tuple val(meta), path('*.interval_list'), emit: interval_list
    path  "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args_file = arguments_file ? "--arguments_file ${arguments_file}" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard BedToIntervalList] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard \\
        -Xmx${avail_mem}g \\
        BedToIntervalList \\
        --INPUT $bed \\
        --OUTPUT ${prefix}.interval_list \\
        --SEQUENCE_DICTIONARY $dict \\
        --TMP_DIR . \\
        $args_file $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard BedToIntervalList --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "picard \\
        -Xmx${avail_mem}g \\
        BedToIntervalList \\
        --INPUT $bed \\
        --OUTPUT ${prefix}.interval_list \\
        --SEQUENCE_DICTIONARY $dict \\
        --TMP_DIR . \\
        $args_file $args"

    touch ${prefix}.interval_list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard BedToIntervalList --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
