// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BWAMEM3_INDEX {
    tag "${meta.id}"
    // NOTE bwa-mem3 builds an FM-index with libsais; peak memory scales with the reference size.
    // `--meth` also builds a converted seed index over twice the text, roughly doubling peak memory.
    memory { 280.MB * Math.ceil(fasta.size() / 10000000) * ("${task.ext.args ?: ''}".contains('--meth') ? 2 : 1) * task.attempt }

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/98/98b0f53f5aed4a3743c17acf6fe59bb0587228a320a6bbe01b3ff3e9ccf85d43/data'
        : 'community.wave.seqera.io/library/bwa-mem3:0.16.0--b56ccc800b165a7e'}"

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
        // `index --meth` rejects -p and writes next to its input, so index a symlink named after the prefix.
        // `mem --meth` needs neither the symlink, the intermediate <prefix>.meth.fa nor the plain <prefix>.bwt.2bit.64
        // (it still reads the plain .amb, .ann and .pac), so a methylation index can't be used for plain alignment.
        """
        mkdir bwamem3
        ln -s ../${fasta} bwamem3/${prefix}
        bwa-mem3 \\
            index \\
            $args \\
            -t ${task.cpus} \\
            bwamem3/${prefix}
        rm bwamem3/${prefix} bwamem3/${prefix}.meth.fa bwamem3/${prefix}.bwt.2bit.64
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
    def index_files = args.contains('--meth') ? "{amb,ann,pac,meth.amb,meth.ann,meth.bwt.2bit.64,meth.pac}" : "{amb,ann,bwt.2bit.64,pac}"
    """
    mkdir bwamem3
    touch bwamem3/${prefix}.${index_files}
    """
}
