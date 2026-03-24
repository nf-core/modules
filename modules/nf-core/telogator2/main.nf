process TELOGATOR2 {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ab/ab4d9d463b2866006f8cbca9fbe6f978b1803e41f2a97d9f4d3c14ff6d97822f/data'
        : 'community.wave.seqera.io/library/telogator2:2.2.3--01b2748e09721f3b' }"

    input:
    tuple val(meta), path(reads), path(reads_index)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("${prefix}"), emit: results
    tuple val("${task.process}"), val('telogator2'), eval("telogator2 --version | sed 's/telogator2 //'"), emit: versions_telogator2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ref_arg = fasta ? "--ref ${fasta}" : ""
    """
    set +e
    telogator2 \\
        -i ${reads} \\
        -o ${prefix} \\
        -p ${task.cpus} \\
        ${ref_arg} \\
        ${args} \\
    > telogator2.log 2>&1
    exit_code=\$?
    cat telogator2.log
    set -e

    if [ \$exit_code -ne 0 ]; then
        if grep -q 'No telomere reads found' telogator2.log; then
            echo "Warning: No telomere reads found in input"
        else
            exit \$exit_code
        fi
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/telomere_length_summary.tsv
    """
}
