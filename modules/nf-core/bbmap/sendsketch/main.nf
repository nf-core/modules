process BBMAP_SENDSKETCH {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bbmap=39.01"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0':
        'biocontainers/bbmap:39.01--h5c4e2a8_0' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("*.txt")  , emit: hits
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_used = file.size() > 1 ? file[0] : file
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[bbmap sendsketch.sh] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    sendsketch.sh -Xmx${avail_mem}M \\
        in=${file_used} \\
        out=${prefix}.txt \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
