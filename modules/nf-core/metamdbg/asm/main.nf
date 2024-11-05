process METAMDBG_ASM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metamdbg:1.0--hdcf5f25_1':
        'biocontainers/metamdbg:1.0--hdcf5f25_1' }"

    input:
    tuple val(meta), path(reads)
    val(input_type)

    output:
    tuple val(meta), path("*.contigs.fasta.gz"), emit: contigs
    tuple val(meta), path("*.metaMDBG.log")    , emit: log
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    switch(input_type) {
        case "hifi": input = "--in-hifi ${reads}"; break
        case "ont" : input = "--in-ont ${reads}" ; break
        default:
            error("ERROR: input_type must be one of either 'hifi' or 'ont'.")
            break
    }
    """
    metaMDBG asm \\
        --threads ${task.cpus} \\
        --out-dir . \\
        ${args} \\
        ${input}

    rm -r tmp/

    mv contigs.fasta.gz ${prefix}.contigs.fasta.gz
    mv metaMDBG.log ${prefix}.metaMDBG.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metamdbg: \$(metaMDBG | grep "Version" | sed 's/ Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.metaMDBG.log
    touch ${prefix}.contigs.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metamdbg: \$(metaMDBG | grep "Version" | sed 's/ Version: //')
    END_VERSIONS
    """
}
