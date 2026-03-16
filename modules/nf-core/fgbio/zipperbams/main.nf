process FGBIO_ZIPPERBAMS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe9479adc5e6e0a1c125d346fdfa0dd313834249e9c55c40e8d44ec3a48c6559/data' :
        'community.wave.seqera.io/library/fgbio:3.1.1--6c9a88faf1d62b6c' }"

    input:
    tuple val(meta), path(unmapped_bam)
    tuple val(meta2), path(mapped_bam)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(dict)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val("${task.process}"), val('fgbio'), eval('fgbio --version 2>&1 | tr -d "[:cntrl:]" | sed -e "s/^.*Version: //;s/\\[.*$//"'), topic: versions, emit: versions_fgbio

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''  // fgbio common options
    def args2  = task.ext.args2 ?: '' // fgbio tool options
    prefix = task.ext.prefix ?: "${meta.id}_zipped"
    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio ZipperBams] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else if (mem_gb > task.memory.giga) {
        if (task.memory.giga < 2) {
            mem_gb = 1
        } else {
            mem_gb = task.memory.giga - 1
        }
    }

    if ("${unmapped_bam}" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    if ("${mapped_bam}" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    if (!args.contains('--async-io=')) {
        args = "--async-io=true ${args}"
    }
    if (!args.contains('--compression ')) {
        args = "--compression 0 ${args}"
    }
    """
    fgbio -Xmx${mem_gb}g \\
        ${args} \\
        ZipperBams \\
        --unmapped ${unmapped_bam} \\
        --input ${mapped_bam} \\
        --ref ${fasta} \\
        ${args2} \\
        --output ${prefix}.bam
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_zipped"
    if ("${unmapped_bam}" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    if ("${mapped_bam}" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    touch ${prefix}.bam
    """
}
