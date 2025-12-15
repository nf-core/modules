process RIBOCODE_PREPARE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribocode:1.2.15--pyhfa5458b_0':
        'biocontainers/ribocode:1.2.15--pyhfa5458b_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("annotation")                                             , emit: annotation
    tuple val("${task.process}"), val('ribocode'), eval('RiboCode --version  2>&1') , emit: versions_ribocode, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    prepare_transcripts \\
        -g ${gtf} \\
        -f ${fasta} \\
        -o annotation \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''

    """
    mkdir annotation

    touch annotation/transcripts_cds.txt
    touch annotation/transcripts_sequence.fa
    touch annotation/transcripts.pickle
    """
}
