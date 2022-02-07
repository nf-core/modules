process PICARD_FILTERSAMREADS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

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
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard FilterSamReads] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    if ( filter == 'includeAligned' || filter == 'excludeAligned' ) {
        """
        picard \\
            FilterSamReads \\
            -Xmx${avail_mem}g \\
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
            -Xmx${avail_mem}g \\
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
}
