// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CONCOCT_EXTRACTFASTABINS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b414706a0ce09cd26055f920ea4d7db8e6e0a4845f806f0918e9792a2c847f95/data'
        : 'community.wave.seqera.io/library/concoct_setuptools:7d62fa8e6704b1fe'}"

    input:
    tuple val(meta), path(original_fasta), path(csv)

    output:
    tuple val(meta), path("${prefix}/*.fa.gz"), emit: fasta
    tuple val("${task.process}"), val('concoct'), eval("concoct --version 2>&1 | sed -n 's/concoct //p'"), topic: versions, emit: versions_concoct

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}

    extract_fasta_bins.py \\
        ${args} \\
        ${original_fasta} \\
        ${csv} \\
        --output_path ${prefix}

    ## Add prefix to each file to disambiguate one sample's 1.fa, 2.fa from sample2
    for i in ${prefix}/*.fa; do
        mv \${i} \${i/\\///${prefix}_}
        gzip \${i/\\///${prefix}_}
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    echo "" | gzip > ${prefix}/${prefix}.fa.gz
    """
}
