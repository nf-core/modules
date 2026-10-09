// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process MISOPY_INDEX {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0a/0a6a585c9b85aee50b8da33d7d4b27d209174ea7c7279a2ac73ceda6b8d4d193/data'
        : 'community.wave.seqera.io/library/python_misopy:9d5a390611c447f5'}"

    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path("${prefix}"), emit: miso_index
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //g"'), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('misopy'), eval('python -c "import pkg_resources; print(pkg_resources.get_distribution(\'misopy\').version)"'), topic: versions, emit: versions_misopy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    grep -v "^#" ${gff3} | \\
    awk 'BEGIN{OFS="\\t"} NF==9 && \$6~/^[0-9.]+\$/ || \$6=="." {print}' > filtered.gff3

    index_gff \\
        --index filtered.gff3 \\
        ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/compressed_ids_to_genes.shelve.{bak,dat,dir}
    touch ${prefix}/genes_to_filenames.shelve.{bak,dat,dir}
    touch ${prefix}/genes.gff
    """
}
