process BBMAP_INDEX {
    tag "$fasta"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
    path fasta

    output:
    path 'ref'                    , emit: index
    tuple val("${task.process}"), val('bbmap'), eval('bbversion.sh | grep -v "Duplicate cpuset"'), emit: versions_bbmap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bbmap.sh \\
        ref=${fasta} \\
        $args \\
        threads=$task.cpus \\
        -Xmx${task.memory.toGiga()}g
    """
    stub:
    """
    mkdir -p ref
    touch ref/info.txt
    touch ref/summary.txt
    """
}
