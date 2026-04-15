process BOWTIE_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6f5ca09fd5aab931d9b87c532c69e0122ce5ff8ec88732f906e12108d48425e9/data' :
        'community.wave.seqera.io/library/bowtie_htslib_samtools:e1e242368ffcb5d3' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    val (save_unaligned)

    output:
    tuple val(meta), path('*.bam')     , emit: bam
    tuple val(meta), path('*.out')     , emit: log
    tuple val(meta), path('*fastq.gz') , emit: fastq, optional : true
    tuple val("${task.process}"), val('bowtie'), eval("bowtie --version 2>&1 | sed -n 's/.*bowtie-align-s version //p'"), emit: versions_bowtie, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions
    tuple val("${task.process}"), val('gzip'), eval("gzip --version 2>&1 | sed '1!d;s/gzip //'"), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def unaligned = save_unaligned ? "--un ${prefix}.unmapped.fastq" : ''
    def endedness = meta.single_end ? "$reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    INDEX=\$(find -L ./ -name "*.3.ebwt" | sed 's/\\.3.ebwt\$//')

    bowtie \\
        --threads $task.cpus \\
        --sam \\
        -x \$INDEX \\
        -q \\
        $unaligned \\
        $args \\
        $endedness \\
        2>| >(tee ${prefix}.out >&2) \\
        | samtools view $args2 -@ $task.cpus -bS -o ${prefix}.bam -

    if [ -f ${prefix}.unmapped.fastq ]; then
        gzip ${prefix}.unmapped.fastq
    fi
    if [ -f ${prefix}.unmapped_1.fastq ]; then
        gzip ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def unaligned = save_unaligned ?
        meta.single_end ? "echo '' | gzip > ${prefix}.unmapped.fastq.gz" :
            "echo '' | gzip > ${prefix}.unmapped_1.fastq.gz; echo '' | gzip > ${prefix}.unmapped_2.fastq.gz" :
                ''
    """
    touch ${prefix}.bam
    touch ${prefix}.out
    $unaligned
    """

}
