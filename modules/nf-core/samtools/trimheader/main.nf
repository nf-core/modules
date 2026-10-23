// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SAMTOOLS_TRIMHEADER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a9/a95118435ca30f8b667ca92f1ba45a26abb58c1f09e96a9e7e8431c73ff9a299/data'
        : 'community.wave.seqera.io/library/htslib_samtools_gawk:6f88a5f6186c9e3d'}"

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
