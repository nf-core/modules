process PARAPHASE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // Needs mulled container with minimap2
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-058de387f9917a7a63953f496cdd203bca83b790:86215829f86df9201683956877a19d025261ff66-0':
        'biocontainers/mulled-v2-058de387f9917a7a63953f496cdd203bca83b790:86215829f86df9201683956877a19d025261ff66-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(config)

    output:
    tuple val(meta), path("*.paraphase.json")              , emit: json
    tuple val(meta), path("*.paraphase.bam")               , emit: bam
    tuple val(meta), path("*.paraphase.bam.bai")           , emit: bai
    tuple val(meta), path("${prefix}_paraphase_vcfs/*.vcf"), emit: vcf, optional: true
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def config_file = config ? "--config $config" : ""
    """
    paraphase \\
        $args \\
        --threads $task.cpus \\
        --bam $bam \\
        --reference $fasta \\
        --prefix $prefix \\
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
    mkdir ${prefix}_paraphase_vcfs

    touch ${prefix}.paraphase.json
    touch ${prefix}.paraphase.bam
    touch ${prefix}.paraphase.bam.bai
    touch ${prefix}_paraphase_vcfs/${prefix}_stub.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paraphase: \$(paraphase --version)
    END_VERSIONS
    """
}
