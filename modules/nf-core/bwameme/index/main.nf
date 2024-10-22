process BWAMEME_INDEX {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa-meme:1.0.6--hdcf5f25_2':
        'biocontainers/bwa-meme:1.0.6--hdcf5f25_2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwameme"), emit: index
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta}"
    def VERSION = '1.0.6' // WARN: Version information provided by tool on CLI is incorrect. Please update this string when bumping container versions.
    """
    mkdir bwameme

    bwa-meme index \\
        $args \\
        -t $task.cpus \\
        -p bwameme/$prefix \\
        $fasta

    build_rmis_dna.sh bwameme/$prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwameme: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta}"
    def VERSION = '1.0.6' // WARN: Version information provided by tool on CLI is incorrect. Please update this string when bumping container versions.
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwameme: $VERSION
    END_VERSIONS
    """
}
