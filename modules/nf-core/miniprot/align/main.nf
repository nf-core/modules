process MINIPROT_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/miniprot:0.11--he4a0461_2':
        'biocontainers/miniprot:0.11--he4a0461_2' }"

    input:
    tuple val(meta), path(pep)
    tuple val(meta2), path(ref)

    output:
    tuple val(meta), path("*.paf"), optional: true, emit: paf
    tuple val(meta), path("*.gff"), optional: true, emit: gff
    tuple val("${task.process}"), val('miniprot'), eval("miniprot --version"), topic: versions, emit: versions_miniprot

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--gff") ? "gff" : "paf"
    """
    miniprot \\
        $args \\
        -t $task.cpus \\
        ${ref} \\
        ${pep} \\
        > ${prefix}.${extension}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--gff") ? "gff" : "paf"
    """
    touch ${prefix}.${extension}
    """
}
