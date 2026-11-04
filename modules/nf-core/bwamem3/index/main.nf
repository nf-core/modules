// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BWAMEM3_INDEX {
    tag "${meta.id}"
    // NOTE bwa-mem3 builds an FM-index with libsais; peak memory scales with the reference size.
    // `--meth` also builds a C->T/G->A converted seed index over twice the text, roughly doubling peak memory.
    memory { 280.MB * Math.ceil(fasta.size() / 10000000) * ("${task.ext.args ?: ''}".contains('--meth') ? 2 : 1) * task.attempt }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/70/70c94630fc1da5d691d57ef22c9c02288838393dadffb1dbf008147441fd77ef/data'
:         'community.wave.seqera.io/library/bwa-mem3:0.15.0--832697eee00ba922' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwamem3"), emit: index
    tuple val("${task.process}"), val('bwamem3'), eval("bwa-mem3 version | sed -nE '1 s/^([0-9]+(\\.[0-9]+)+).*/\\1/p'"), emit: versions_bwamem3, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta}"
    def args = task.ext.args ?: ''
    if (args.contains('--meth')) {
        // `index --meth` rejects -p and writes <in.fasta>.* and <in.fasta>.meth.* next to its input,
        // so index a symlink inside the output directory. Drop the symlink and the intermediate
        // converted FASTA (<in.fasta>.meth.fa) afterwards: `mem --meth` needs neither.
        """
        mkdir bwamem3
        ln -s ../${fasta} bwamem3/${prefix}
        bwa-mem3 \\
            index \\
            $args \\
            -t ${task.cpus} \\
            bwamem3/${prefix}
        rm bwamem3/${prefix} bwamem3/${prefix}.meth.fa
        """
    } else {
        """
        mkdir bwamem3
        bwa-mem3 \\
            index \\
            $args \\
            -t ${task.cpus} \\
            -p bwamem3/${prefix} \\
            $fasta
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${fasta}"
    def args = task.ext.args ?: ''
    def meth_files = args.contains('--meth') ? "touch bwamem3/${prefix}.meth.{amb,ann,bwt.2bit.64,pac}" : ''
    """
    mkdir bwamem3
    touch bwamem3/${prefix}.0123
    touch bwamem3/${prefix}.amb
    touch bwamem3/${prefix}.ann
    touch bwamem3/${prefix}.bwt.2bit.64
    touch bwamem3/${prefix}.pac
    ${meth_files}
    """
}
