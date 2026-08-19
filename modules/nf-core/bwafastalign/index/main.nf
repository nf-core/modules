process BWAFASTALIGN_INDEX {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f8/f8c975324a12014c8a817c2c1ad0cd68b077cf09c4370717589177b262dcd1dc/data':
        'community.wave.seqera.io/library/bwa-fastalign_mbuffer_samtools:35f24ce8addcd26b'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwafastalign"), emit: index
    tuple val("${task.process}"), val('bwafastalign'), val('1.0.0'), topic: versions, emit: versions_bwafastalign
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    tuple val("${task.process}"), val('samtools'), eval("samtools --version 2>&1 | sed '1!d;s/.* //'") , topic: versions, emit: versions_samtools
    tuple val("${task.process}"), val('mbuffer'), eval("mbuffer --version 2>&1 | sed -n 's/mbuffer version //p'") , topic: versions, emit: versions_mbuffer

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta}"
    """
    mkdir bwafastalign

    bwa-fastalign index \\
        $args \\
        -p bwafastalign/$prefix \\
        $fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta}"
    """
    mkdir bwafastalign
    touch bwafastalign/${prefix}.amb
    touch bwafastalign/${prefix}.ann
    touch bwafastalign/${prefix}.bwt
    touch bwafastalign/${prefix}.bytesa
    touch bwafastalign/${prefix}.fmt
    touch bwafastalign/${prefix}.kmer
    touch bwafastalign/${prefix}.pac
    touch bwafastalign/${prefix}.sa
    """
}
