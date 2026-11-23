// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BISMARK_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/20/20bce4b3bce2e7045431676190322012e22a96cef1a53eaadd14acc89de2e232/data'
        : 'community.wave.seqera.io/library/bismark_bwa-mem3:1eb534160cb7def7'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta, stageAs: 'tmp/*')
    // This change mounts as directory containing the FASTA file to prevent nested symlinks
    tuple val(meta3), path(index)
    tuple val(meta4), path(bwamem3_index, stageAs: 'bwamem3_index/*')

    output:
    tuple val(meta), path("*bam"), emit: bam
    tuple val(meta), path("*report.txt"), emit: report
    tuple val(meta), path("*fq.gz"), emit: unmapped, optional: true
    tuple val("${task.process}"), val("bismark"), eval("bismark --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+'"), topic: versions, emit: versions_bismark
    tuple val("${task.process}"), val('bwamem3'), eval("bwa-mem3 version | sed -nE '1 s/^([0-9]+(\\.[0-9]+)+).*/\\1/p'"), emit: versions_bwamem3, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if (task.ext.prefix) {
        args += " --prefix ${task.ext.prefix}"
    }
    def fastq = meta.single_end ? reads : "-1 ${reads[0]} -2 ${reads[1]}"

    // 5-Base and TAPS align to the unconverted genome, so the staged FASTA folder can stand in for the index
    def bwamem3_mode = (args =~ /(^|\s)--(bwamem3|taps)(\s|$)/).find()
    def five_base = bwamem3_mode || (args =~ /(^|\s)--(illumina_5base|five_base)(\s|$)/).find()
    if (bwamem3_index && !bwamem3_mode) {
        error("bismark/align: bwamem3_index needs --illumina_5base --bwamem3 or --taps in ext.args")
    }
    if (bwamem3_index && args.contains('--five_base_index')) {
        error("bismark/align: give the bwa-mem3 index as the bwamem3_index input or as --five_base_index in ext.args, not both")
    }
    if (bwamem3_mode && !index && !bwamem3_index && !args.contains('--five_base_index')) {
        error("bismark/align: --bwamem3/--taps needs bwamem3_index, --five_base_index in ext.args, or an index folder holding Bwamem3_Index/")
    }
    def genome = index ?: (five_base ? 'tmp' : '')
    if (!genome) {
        error("bismark/align: index is required")
    }
    def resolve_index = bwamem3_index ? """
    PREFIX=\$(find -L bwamem3_index/ -name '*.meth.amb' | sed 's/\\.meth\\.amb\$//')
    if [ "\$(printf '%s\\n' "\$PREFIX" | grep -c .)" -ne 1 ]; then
        echo "bwamem3_index must hold exactly one bwa-mem3 --meth index (*.meth.amb), found: \${PREFIX:-none}" >&2
        exit 1
    fi
    """ : ''
    def index_arg = bwamem3_index ? '--five_base_index "$PREFIX"' : ''

    // Give each strand instance `-p = cpus / n_instances` threads (one index load per strand,
    // byte-identical) instead of forking N --multicore chunks that each re-load the index.
    // minimap2/rammap can't take -p from bismark yet, so they keep the legacy --multicore path.
    def isMinimapLike = args.contains('--minimap2') || args.contains('--mm2') || args.contains('--rammap') || args.contains('--ram')
    if (isMinimapLike) {

        // Legacy --multicore auto-compute (unchanged) — minimap2/rammap only.
        if (!args.contains('--multicore') && task.cpus) {

            // Numbers based on recommendation by Felix for a typical mouse genome
            def ccore = 1
            def cpu_per_multicore = 3
            def mem_per_multicore = 13.GB.toBytes()
            if (args.contains('--non_directional')) {
                cpu_per_multicore = 5
                mem_per_multicore = 18.GB.toBytes()
            }

            // How many multicore splits can we afford with the cpus we have?
            ccore = ((task.cpus as int) / cpu_per_multicore) as int

            // Check that we have enough memory
            if (task.memory) {
                def mcore = (task.memory.toBytes() / mem_per_multicore) as int
                ccore = Math.min(ccore, mcore)
            }
            else {
                log.warn("Not able to define bismark align multicore based on available memory")
            }
            if (ccore > 1) {
                args += " --multicore ${ccore}"
            }
        }
    }
    else if (bwamem3_mode) {
        // One bwa-mem3 instance; bismark requires -p >= 2
        if (!(args =~ /(^|\s)(--multicore|--parallel|--threads|-p)(\s|=|\d|$)/).find()) {
            args += " -p ${Math.max(2, (task.cpus ?: 2) as int)}"
        }
    }
    else if (!args.contains('--multicore') && !args.contains('--parallel') && !args.contains('-p ') && task.cpus) {

        def n_instances = 2
        if (args.contains('--combined_index')) {
            // combined index = one both-strands pass; parallel non-directional runs two.
            n_instances = args.contains('--non_directional') && args.contains('--combined_index_parallel') ? 2 : 1
        }
        else if (args.contains('--non_directional')) {
            n_instances = 4
        }

        // bismark requires -p >= 2, so omit it below that (one thread per instance).
        def pthreads = ((task.cpus as int) / n_instances) as int
        if (pthreads >= 2) {
            args += " -p ${pthreads}"
        }
    }
    """
    ${resolve_index}
    bismark \\
        ${fastq} \\
        --genome ${genome} \\
        --bam \\
        ${index_arg} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.report.txt
    """
}
