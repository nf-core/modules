process SAMTOOLS_COLLATEFASTQ {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(fasta), path(fai)
    val interleave

    output:
    tuple val(meta), path("*_{1,2}.fq.gz"), emit: fastq, optional: true
    tuple val(meta), path("*_interleaved.fq"), emit: fastq_interleaved, optional: true
    tuple val(meta), path("*_other.fq.gz"), emit: fastq_other
    tuple val(meta), path("*_singleton.fq.gz"), emit: fastq_singleton, optional: true
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    def output = interleave && !meta.single_end
        ? "> ${prefix}_interleaved.fq"
        : meta.single_end
            ? "-1 ${prefix}_1.fq.gz -s ${prefix}_singleton.fq.gz"
            : "-1 ${prefix}_1.fq.gz -2 ${prefix}_2.fq.gz -s ${prefix}_singleton.fq.gz"

    """
    samtools collate \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reference} \\
        -O \\
        ${input} \\
        . |

    samtools fastq \\
        ${args2} \\
        --threads ${task.cpus} \\
        ${reference} \\
        -0 ${prefix}_other.fq.gz \\
        ${output}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def empty = "echo '' | gzip "
    def singletoncommand = "${empty}> ${prefix}_singleton.fq.gz"
    def interleavecommand = interleave && !meta.single_end ? "${empty}> ${prefix}_interleaved.fq.gz" : ""
    def output1command = !interleave ? "${empty}> ${prefix}_1.fq.gz" : ""
    def output2command = !interleave && !meta.single_end ? "${empty}> ${prefix}_2.fq.gz" : ""

    """
    ${output1command}
    ${output2command}
    ${interleavecommand}
    ${singletoncommand}
    ${empty}> ${prefix}_other.fq.gz
    """
}
