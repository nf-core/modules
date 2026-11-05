// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BWAMETH_INDEX {
    tag "${fasta}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/44/4405ebf9fd0b415a9e28b0cb503b11af0e86c88ab335b098ff4cf48f8f0162e1/data'
:         'community.wave.seqera.io/library/bwameth:0.2.10--a589797f382631bc' }"

    input:
    tuple val(meta), path(fasta, name: "BwamethIndex/")
    val use_mem2

    output:
    tuple val(meta), path("BwamethIndex"), emit: index
    tuple val("${task.process}"), val('bwameth'), eval("bwameth.py --version | cut -f2 -d' '"), emit: versions_bwameth_index, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def index_cmd = use_mem2 ? "index-mem2" : "index"
    """
    bwameth.py ${index_cmd} ${fasta}

    rm $fasta
    """

    stub:
    """
    rm ${fasta}

    mkdir -p BwamethIndex/
    touch BwamethIndex/genome.fasta.bwameth.c2t
    touch BwamethIndex/genome.fasta.bwameth.c2t.amb
    touch BwamethIndex/genome.fasta.bwameth.c2t.ann
    touch BwamethIndex/genome.fasta.bwameth.c2t.bwt
    touch BwamethIndex/genome.fasta.bwameth.c2t.pac
    touch BwamethIndex/genome.fasta.bwameth.c2t.sa
    """
}
