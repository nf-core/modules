process PARAPHASE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // Needs mulled container with minimap2
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-058de387f9917a7a63953f496cdd203bca83b790:86215829f86df9201683956877a19d025261ff66-0':
        'biocontainers/mulled-v2-058de387f9917a7a63953f496cdd203bca83b790:86215829f86df9201683956877a19d025261ff66-0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.json")              , emit: json
    tuple val(meta), path("*.bam")               , emit: bam
    tuple val(meta), path("*.bai")               , emit: bai
    tuple val(meta), path("${prefix}_vcfs/*.vcf"), emit: vcf
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    paraphase \\
        $args \\
        --threads $task.cpus \\
        --bam $bam \\
        --reference $fasta \\
        --prefix $prefix \\
        --gene DDT \\
        --out .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paraphase: \$(paraphase --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}_vcfs

    touch ${prefix}.json
    touch ${prefix}_realinged_tagged.bam
    touch ${prefix}_realigned_tagged.bam.bai
    touch ${prefix}_vcfs/${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paraphase: \$(paraphase --version)
    END_VERSIONS
    """
}
