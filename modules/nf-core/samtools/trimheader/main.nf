// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SAMTOOLS_TRIMHEADER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bf/bf776bd8977ad787661a64512dbecfbe70f88b504f939b473aa4ad87b1747edd/data'
        : 'community.wave.seqera.io/library/htslib_samtools_gawk:5366169842fe8538'}"

    input:
    tuple val(meta), path(bam, stageAs: 'input/*'), path(bai, stageAs: 'input/*')

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools idxstats --threads ${task.cpus} ${bam} > idxstats.tsv
    {
        samtools view -H ${bam} \\
            | awk -F '\\t' 'NR == FNR { if (\$3 + \$4 > 0) keep[\$1]; next } !/^@SQ/ || (substr(\$2, 4) in keep)' idxstats.tsv -
        samtools view --threads ${task.cpus} ${bam}
    } | samtools view ${args} --threads ${task.cpus} -b -o ${prefix}.bam -
    rm idxstats.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
