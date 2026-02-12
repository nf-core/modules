process DEACON_INDEX {
    tag "$fasta"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deacon:0.13.2--h7ef3eeb_1':
        'biocontainers/deacon:0.13.2--h7ef3eeb_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.idx"), emit: index
    tuple val("${task.process}"), val('deacon'), eval('deacon --version | head -n1 | sed "s/deacon //g"'), emit: versions_deacon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    deacon \\
        index \\
        build \\
        --threads ${task.cpus} \\
        $args \\
        $fasta > ${prefix}.idx
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    touch ${prefix}.idx
    """
}
