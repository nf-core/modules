
process WHATSHAP_HAPLOTAG {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d837709891c2d98fc0956f6fd0dba18b0f67d96c4db74ccbae7db98fd00afe42/data' :
        'community.wave.seqera.io/library/whatshap:2.8--7fe530bc624a3e5a' }"

    input:
    tuple val(meta),  path(vcf), path(tbi), path(bam), path(bai)
    tuple val(meta2), path(fasta) // empty channel [] if not needed
    tuple val(meta3), path(fai) // empty channel [] if not needed
    val(include_tsv_output)    // value:   [ true | false ]

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.{tsv,tsv.gz}"), emit: tsv, optional: true
    tuple val("${task.process}"), val('whatshap'), eval("whatshap --version"), topic: versions, emit: versions_whatshap

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_haplotagged"

    def reference = fasta ? "--reference $fasta" : "--no-reference"
    def output_tsv = include_tsv_output  ? "--output-haplotag-list ${prefix}.tsv.gz" : ''

    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    whatshap haplotag \\
        $args \\
        -o ${prefix}.bam \\
        $reference \\
        --output-threads $task.cpus \\
        $output_tsv \\
        ${vcf} \\
        ${bam}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def output_tsv = include_tsv_output  ? "echo '' | gzip > ${prefix}.tsv.gz" : ''
    """
    touch ${prefix}.bam
    $output_tsv

    echo $args
    """
}
