process PICARD_MERGESAMFILES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::picard=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1' :
        'quay.io/biocontainers/picard:3.0.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_files = bams.sort()
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard MergeSamFiles] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    if (bam_files.size() > 1) {
        """
        picard \\
            -Xmx${avail_mem}M \\
            MergeSamFiles \\
            $args \\
            ${'--INPUT '+bam_files.join(' --INPUT ')} \\
            --OUTPUT ${prefix}.bam
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            picard: \$( echo \$(picard MergeSamFiles --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
        END_VERSIONS
        """
    } else {
        """
        ln -s ${bam_files[0]} ${prefix}.bam
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            picard: \$( echo \$(picard MergeSamFiles --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
        END_VERSIONS
        """
    }
}
