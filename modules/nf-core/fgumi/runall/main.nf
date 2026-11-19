// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FGUMI_RUNALL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bc/bc7499ff550a8302bcb8c75d35275b45519dc42f5e810fb078c50db6ce5c54e9/data'
        : 'community.wave.seqera.io/library/bwa-mem3_fgumi_findutils:334cfaafb9102360'}"

    input:
    tuple val(meta), path(reads), val(library)
    tuple val(meta2), path(fasta), path(fai), path(dict)
    tuple val(meta3), path(index)
    val start_from
    val stop_after

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi
    tuple val("${task.process}"), val('bwamem3'), eval("bwa-mem3 version | sed -nE '1 s/^([0-9]+(\\.[0-9]+)+).*/\\1/p'"), topic: versions, emit: versions_bwamem3

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input = start_from == 'extract'
        ? "--extract::inputs ${reads} --extract::sample ${meta.id} --extract::library ${library}"
        : "--input ${reads}"
    if ("${reads}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    // The aligner reads its index from the `--ref` path, so the FASTA and the index must share one prefix
    def link_reference = index
        ? """
        INDEX=`find -L ${index}/ -name "*.amb" ! -name "*.meth.amb" | sed 's/\\.amb\$//'`
        REF=reference/`basename \$INDEX`
        mkdir reference
        for file in \$INDEX.*; do ln -s ../\$file reference/; done
        ln -s ../${fasta} \$REF
        ln -s ../${fai} \$REF.fai
        ln -s ../${dict} \$REF.dict
        """.stripIndent()
        : "REF=${fasta}"
    def reference = fasta ? "--ref \$REF" : ''
    """
    ${link_reference}

    fgumi runall \\
        --start-from ${start_from} \\
        --stop-after ${stop_after} \\
        ${input} \\
        ${reference} \\
        --output ${prefix}.bam \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
