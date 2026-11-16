// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process QUALIMAP_BAMQCCRAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a67cd2fed2d1e951052fdea57c914851c08a3d8e32a11937636d1d84ba36b201/data'
:         'community.wave.seqera.io/library/qualimap_samtools:36c0048b499c176d' }"

    input:
    tuple val(meta), path(cram), path(crai)
    path  gff
    path  fasta
    path  fasta_fai

    output:
    tuple val(meta), path("${prefix}"), emit: results
    tuple val("${task.process}"), val('qualimap'), eval("qualimap -h | sed -n 's/^QualiMap v.//p'"), topic: versions, emit: versions_qualimap
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    def collect_pairs = meta.single_end ? '' : '--collect-overlap-pairs'
    def memory = (task.memory.mega*0.8).intValue() + 'M'
    def regions = gff ? "--gff ${gff}" : ''

    """
    unset DISPLAY
    mkdir -p tmp
    export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp

    samtools view -hb -T ${fasta} ${cram} |
    qualimap \\
        --java-mem-size=${memory} \\
        bamqc \\
        ${args} \\
        -bam /dev/stdin \\
        ${regions} \\
        ${collect_pairs} \\
        -outdir ${prefix} \\
        -nt ${task.cpus}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    """
}
