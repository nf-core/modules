process BWAMEME_INDEX {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9d/9ddd41b93c5e182db9d643ca266dd1677e59593a9cb49904b982ff45ad5aa8c3/data':
        'community.wave.seqera.io/library/bwa-meme_mbuffer_samtools:03f3f60b6c289776' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwameme"), emit: index
    tuple val("${task.process}"), val('bwameme'), val('1.0.6'), topic: versions, emit: versions_bwameme
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta}"
    """
    mkdir bwameme

    bwa-meme index \\
        $args \\
        -t $task.cpus \\
        -p bwameme/$prefix \\
        $fasta

    build_rmis_dna.sh bwameme/$prefix
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta}"
    """
    mkdir bwameme
    touch bwameme/${prefix}.0123
    touch bwameme/${prefix}.ann
    touch bwameme/${prefix}.pac
    touch bwameme/${prefix}.amb
    touch bwameme/${prefix}.pos_packed
    touch bwameme/${prefix}.suffixarray_uint64
    touch bwameme/${prefix}.suffixarray_uint64_L0_PARAMETERS
    touch bwameme/${prefix}.suffixarray_uint64_L1_PARAMETERS
    touch bwameme/${prefix}.suffixarray_uint64_L2_PARAMETERS
    """
}
