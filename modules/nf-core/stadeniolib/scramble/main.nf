// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process STADENIOLIB_SCRAMBLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/de/deffb1131242d30a8e6f760a5f9ee85846b0cc74563bcea3fd22c1743ef396fe/data'
        : 'community.wave.seqera.io/library/staden_io_lib:1.15.1--cdbf259bf64228b4'}"

    input:
    tuple val(meta), path(reads)
    path(fasta)
    path(fai)
    path(gzi)

    output:
    tuple val(meta), path("*.{cram,bam}")   , emit: cram
    path "*.gzi"                            , emit: gzi         , optional: true
    tuple val("${task.process}"), val('scramble'), eval("scramble -h | head -n 1 | sed 's/^.*version //'"), emit: versions_scramble, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def inputformat = reads.extension
    def outputformat = "cram"
    if ( "-O sam" in args ) {
        outputformat = "sam"
    } else if ( "-O bam" in args ) {
        outputformat = "bam"
    }

    def reference = ( fasta && fai ) ? "-r ${fasta}" : ''
    if ( outputformat == "cram" && !reference ) {
        error "Cannot convert to CRAM without a reference"
    }

    def gz_index = gzi ? "-g ${gzi}" : ''

    if ( ( outputformat == "cram" || outputformat == "sam" ) && gz_index ) {
        gz_index = ''
        error "Cannot use gzip index for CRAM or SAM output"
    }

    """
    scramble \
        $args \
        -I ${inputformat} \
        $reference \
        -t $task.cpus \
        ${reads} \
        ${prefix}.${outputformat}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def outputformat = "cram"
    if ( "-O sam" in args ) {
        outputformat = "sam"
    } else if ( "-O bam" in args ) {
        outputformat = "bam"
    }

    def reference = ( fasta && fai ) ? "-r ${fasta}" : ''
    if ( outputformat == "cram" && !reference ) {
        error "Cannot convert to CRAM without a reference"
    }

    def gz_index = gzi ? "-g ${gzi}" : ''

    if ( ( outputformat == "cram" || outputformat == "sam" ) && gz_index ) {
        gz_index = ''
        error "Cannot use gzip index for CRAM or SAM output"
    }

    """
    touch ${prefix}.${outputformat}
    """
}
