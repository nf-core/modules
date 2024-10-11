process OPENMS_IDMASSACCURACY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.2.0--haddbca4_4' :
        'biocontainers/openms:3.2.0--haddbca4_4' }"

    input:
    tuple val(meta), path(mzmls), path(idxmls)

    output:
    tuple val(meta), path("*frag_mass_err.tsv") , emit: frag_err
    tuple val(meta), path("*prec_mass_err.tsv") , emit: prec_err, optional: true
    path "versions.yml"                         , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        OpenMS: \$(FileInfo 2>&1 | grep -E '^Version(.*)' | cut -d ' ' -f 2 | cut -d '-' -f 1)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_frag_mass_err.tsv
    touch ${prefix}_prec_mass_err.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        OpenMS: \$(FileInfo 2>&1 | grep -E '^Version(.*)' | cut -d ' ' -f 2 | cut -d '-' -f 1)
    END_VERSIONS
    """
}
