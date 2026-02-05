process METAMDBG_ASM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metamdbg:1.2--h077b44d_0':
        'biocontainers/metamdbg:1.2--h077b44d_0' }"

    input:
    tuple val(meta), path(reads, arity: '1..*')
    val(input_type)

    output:
    tuple val(meta), path("*.contigs.fasta.gz"), emit: contigs
    tuple val(meta), path("*.metaMDBG.log")    , emit: log
    tuple val("${task.process}"), val('metamdbg'), eval('metaMDBG | sed -n "s/.*Version: //p"'), emit: versions_metamdbg, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if(!(input_type in ["hifi", "ont"])) {
        error("ERROR: input_type must be one of either 'hifi' or 'ont'.")
    }
    """
    metaMDBG asm \\
        --threads ${task.cpus} \\
        --out-dir . \\
        ${args} \\
        --in-${input_type} ${reads}

    rm -r tmp/

    mv contigs.fasta.gz ${prefix}.contigs.fasta.gz
    mv metaMDBG.log ${prefix}.metaMDBG.log
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    touch ${prefix}.metaMDBG.log
    touch ${prefix}.contigs.fasta.gz
    """
}
