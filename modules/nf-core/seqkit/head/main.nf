// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SEQKIT_HEAD {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/29/297e985dd395d87f9089b314dbe685ddce59fa0770c4d1f5e1009bc7b6ec7a8d/data'
        : 'community.wave.seqera.io/library/seqkit:2.13.0--3bacfb68370dd835'}"

    input:
    tuple val(meta), path(fastqs), val(seq_count)

    output:
    tuple val(meta), path("${prefix}_subset_*"), emit: subset
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), emit: versions_seqkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    for f in ${fastqs.join(' ')}
    do
        seqkit head \\
            ${args} \\
            --threads ${task.cpus} \\
            -n ${seq_count} \\
            -o "${prefix}_subset_\$(basename \$f)" \\
            \$f
    done
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    for f in ${fastqs.join(' ')}
    do
        echo '' | gzip > "${prefix}_subset_\$(basename \$f)"
    done
    """
}
