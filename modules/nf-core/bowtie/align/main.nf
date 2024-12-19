process BOWTIE_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8c0819a9b1f520c49c933e667ae50de2a0730ece4c8b9efe79ac5e403963a9f/data' :
        'community.wave.seqera.io/library/bowtie_samtools:e1a14e1ce4e0170d' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    val (save_unaligned)

    output:
    tuple val(meta), path("*.bam")     , emit: bam
    tuple val(meta), path("*.out")     , emit: log
    tuple val(meta), path("*fastq.gz") , emit: fastq, optional : true
    path  "versions.yml"               , emit: versions

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
        2> >(tee ${prefix}.out >&2) \\
        | samtools view $args2 -@ $task.cpus -bS -o ${prefix}.bam -

    if [ -f ${prefix}.unmapped.fastq ]; then
        gzip ${prefix}.unmapped.fastq
    fi
    if [ -f ${prefix}.unmapped_1.fastq ]; then
        gzip ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def unaligned = save_unaligned ?
                    meta.single_end ? "echo '' | gzip > ${prefix}.unmapped.fastq.gz" :
                        "echo '' | gzip > ${prefix}.unmapped_1.fastq.gz; echo '' | gzip > ${prefix}.unmapped_2.fastq.gz"
                    : ''
    """
    touch ${prefix}.bam
    touch ${prefix}.out
    $unaligned

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """


}
