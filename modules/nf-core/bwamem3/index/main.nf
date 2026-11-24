// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BWAMEM3_INDEX {
    tag "${meta.id}"
    // NOTE bwa-mem3 builds an FM-index with libsais; peak memory scales with the reference size.
    // `--meth` also builds a converted seed index over twice the text, roughly doubling peak memory.
    memory { 280.MB * Math.ceil(fasta.size() / 10000000) * ("${task.ext.args ?: ''}".contains('--meth') ? 2 : 1) * task.attempt }

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d6/d69df8cc3749b0a3a3f687d76d9772beff6ef1aa75e2bcf47c451b6c86bbed88/data'
        : 'community.wave.seqera.io/library/bwa-mem3:0.18.0--a6470642c1454740'}"

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
    // Without --max-memory, bwa-mem3 holds back max(2 GB, 5%) of the memory it detects, up to half,
    // which leaves too little of a small task's request (e.g. a --meth index of human chr21).
    def max_memory = task.memory && !args.contains('--max-memory') ? "--max-memory ${task.memory.toMega()}M" : ''
    // A --meth build skips the plain FM-index, which `mem --meth` never reads, and writes its scratch FASTA to the work dir.
    def meth_only = args.contains('--meth') && !args.contains('--meth-only') ? '--meth-only' : ''
    def tmp_dir = args.contains('--meth') && !args.contains('--tmp-dir') ? '--tmp-dir .' : ''
    """
    mkdir bwamem3
    bwa-mem3 \\
        index \\
        $args \\
        ${meth_only} \\
        ${tmp_dir} \\
        ${max_memory} \\
        -t ${task.cpus} \\
        -p bwamem3/${prefix} \\
        $fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta}"
    def args = task.ext.args ?: ''
    def index_files = args.contains('--meth') ? "{amb,ann,pac,meth.amb,meth.ann,meth.bwt.2bit.64,meth.pac}" : "{amb,ann,bwt.2bit.64,pac}"
    """
    mkdir bwamem3
    touch bwamem3/${prefix}.${index_files}
    """
}
