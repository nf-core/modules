process BBMAP_SENDSKETCH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'wave.seqera.io/wt/5ce851926d93/wave/build:bbmap-39.18_pigz-2.8--4a8254d490b5baa0' }"

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
