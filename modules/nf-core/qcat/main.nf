// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process QCAT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/45/457da82216e5355e298ecd5b4ba3d5c5a81421605af9029da08680c0e31f2ba4/data'
        : 'community.wave.seqera.io/library/qcat_python_setuptools:f1f481513e1ca893'}"

    input:
    tuple val(meta), path(reads)
    val barcode_kit

    output:
    tuple val(meta), path("fastq/*.fastq.gz"), emit: reads
    tuple val("${task.process}"), val('qcat'), eval("qcat --version | sed 's/.*qcat //;s/ .*//'"), emit: versions_qcat, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ## Unzip fastq file
    ## qcat doesn't support zipped files yet
    FILE=${reads}
    if [[ \$FILE == *.gz ]]
    then
        zcat ${reads} > unzipped.fastq
        FILE=unzipped.fastq
    fi

    qcat \\
        -f \$FILE \\
        -b ./fastq \\
        --kit ${barcode_kit}

    ## Zip fastq files
    gzip fastq/*
    """

    stub:
    """
    mkdir -p fastq

    echo "" | gzip > fastq/barcode00.fastq.gz
    echo "" | gzip > fastq/none.fastq.gz
    """
}
