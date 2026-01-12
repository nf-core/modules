process OPENMS_IDMASSACCURACY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.4.1--h81ffffe_1' :
        'biocontainers/openms:3.4.1--h81ffffe_1' }"

    input:
    tuple val(meta), path(mzmls), path(idxmls)

    output:
    tuple val(meta), path("*frag_mass_err.tsv") , emit: frag_err
    tuple val(meta), path("*prec_mass_err.tsv") , emit: prec_err, optional: true
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | grep -E '^Version' | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    IDMassAccuracy \\
        -in $mzmls \\
        -id_in $idxmls \\
        -out_fragment ${prefix}_frag_mass_err.tsv \\
        -threads $task.cpus \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_frag_mass_err.tsv
    touch ${prefix}_prec_mass_err.tsv
    """
}
