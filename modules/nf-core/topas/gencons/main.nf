process TOPAS_GENCONS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::topas=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/topas:1.0.1--hdfd78af_1':
        'biocontainers/topas:1.0.1--hdfd78af_1' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta), path(reference)

    output:
    tuple val(meta), path("*.fasta.gz"), emit: fasta
    tuple val(meta), path("*.ccf")     ,      emit: ccf
    tuple val(meta), path("*.log")     ,      emit: log
    path "versions.yml"                ,      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """

    topas \\
        GenConS \\
        $args \\
        -o ${prefix}.fasta \\
        -snps $vcf \\
        -ref $reference

    gzip ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        topas: $VERSION
    END_VERSIONS
    """
}
