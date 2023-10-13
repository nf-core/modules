process PICARD_FILTERSAMREADS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::picard=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1' :
        'biocontainers/picard:3.0.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(bam), path(readlist)
    val filter

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard FilterSamReads] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    if ( filter == 'includeAligned' || filter == 'excludeAligned' ) {
        """
        picard \\
            FilterSamReads \\
            -Xmx${avail_mem}M \\
            --INPUT $bam \\
            --OUTPUT ${prefix}.bam \\
            --FILTER $filter \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            picard: \$(picard FilterSamReads --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
        END_VERSIONS
        """
    } else if ( filter == 'includeReadList' || filter == 'excludeReadList' ) {
        """
        picard \\
            FilterSamReads \\
            -Xmx${avail_mem}M \\
            --INPUT $bam \\
            --OUTPUT ${prefix}.bam \\
            --FILTER $filter \\
            --READ_LIST_FILE $readlist \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            picard: \$(picard FilterSamReads --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard FilterSamReads --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

}
