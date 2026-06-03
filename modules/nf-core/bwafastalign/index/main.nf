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

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta}"
    """
    mkdir bwafastalign

    bwa-fastalign index \\
        $args \\
        -t $task.cpus \\
        -p bwafastalign/$prefix \\
        $fasta

    build_rmis_dna.sh bwafastalign/$prefix
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta}"
    """
    mkdir bwafastalign
    touch bwafastalign/${prefix}.0123
    touch bwafastalign/${prefix}.ann
    touch bwafastalign/${prefix}.pac
    touch bwafastalign/${prefix}.amb
    touch bwafastalign/${prefix}.pos_packed
    touch bwafastalign/${prefix}.suffixarray_uint64
    touch bwafastalign/${prefix}.suffixarray_uint64_L0_PARAMETERS
    touch bwafastalign/${prefix}.suffixarray_uint64_L1_PARAMETERS
    touch bwafastalign/${prefix}.suffixarray_uint64_L2_PARAMETERS
    """
}
