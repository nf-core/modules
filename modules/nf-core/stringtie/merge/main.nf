process STRINGTIE_MERGE {
    tag "$meta.id"
    label 'process_medium'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'biocontainers/stringtie:2.2.1--hecb563c_2' }"

    input:
    tuple val(meta), path(bam)
    path annotation_gtf

    output:
    tuple val(meta), path("*.merged.gtf"), emit: gtf
    tuple val("${task.process}"), val('stringtie'), eval('stringtie --version'), emit: versions_stringtie, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def reference = annotation_gtf ? "-G $annotation_gtf" : ""
    """
    stringtie \\
        --merge \\
        $bam \\
        $reference \\
        -o ${prefix}.merged.gtf \\
        -p $task.cpus \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.gtf
    """
}
