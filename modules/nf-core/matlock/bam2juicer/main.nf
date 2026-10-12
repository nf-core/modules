process MATLOCK_BAM2JUICER {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update the VERSION variable when bumping
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/59/59df5236cc790cac47380a5654f39052a1cb3f9c7868ed397e4b3205a9fb2776/data':
        'community.wave.seqera.io/library/matlock_samtools:3c30bc2808902fde' }"

    input:
    tuple val(meta), path(input_file)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.links.txt")    , emit: links_txt
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools
    tuple val("${task.process}"), val('matlock'), val('20181227'), topic: versions, emit: versions_matlock
    // WARN: Version information not provided by tool on CLI. Please update the VERSION variable when bumping

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.links"
    def args2 = task.ext.args2 ?: ''
    def input_ext = input_file.extension
    def reference = fasta ? "--reference ${fasta}" : ""
    def filter_cmd = args2 != '' ? "samtools view --threads $task.cpus $args2 -Sb $input_file $reference > ${prefix}.juicer.bam" : ''
    def juicer_input = args2 != '' ? "${prefix}.juicer.${input_ext}" : "$input_file"
    def bam_del_cmd = args2 != '' ? "rm ${prefix}.juicer.bam" : ''
    """
    $filter_cmd

    matlock \\
        bam2 \\
        juicer \\
        $juicer_input \\
        ${prefix}.txt

    $bam_del_cmd
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.links"
    def args2 = task.ext.args2 ?: ''
    def input_ext = input_file.extension
    def reference = fasta ? "--reference ${fasta}" : ""
    def filter_cmd = args2 != '' ? "samtools view --threads $task.cpus $args2 -Sb $input_file $reference > ${prefix}.juicer.bam" : ''
    def juicer_input = args2 != '' ? "${prefix}.juicer.${input_ext}" : "$input_file"
    def bam_del_cmd = args2 != '' ? "rm ${prefix}.juicer.bam" : ''
    """
    echo "filter command: $filter_cmd"
    echo "juicer input: $juicer_input"
    echo "bam delete command: $bam_del_cmd"

    touch ${prefix}.txt
    """
}
