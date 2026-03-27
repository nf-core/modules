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
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
    WARNING: This module has been deprecated. Please use nf-core/modules/path/to/new/module

    Reason:
    This module is no longer fit for purpose because not part of openms 3.5.0 version anymore
    """
    assert false: deprecation_message

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
        OpenMS: \$(FileInfo --help 2>&1 | sed -nE 's/^Version: ([0-9.]+).*/\\1/p')
    END_VERSIONS
    """

    stub:
    def deprecation_message = """
    WARNING: This module has been deprecated. Please use nf-core/modules/path/to/new/module

    Reason:
    This module is no longer fit for purpose because not part of openms 3.5.0 version anymore
    """
    assert false: deprecation_message
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_frag_mass_err.tsv
    touch ${prefix}_prec_mass_err.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        OpenMS: \$(FileInfo --help 2>&1 | sed -nE 's/^Version: ([0-9.]+).*/\\1/p')
    END_VERSIONS
    """
}
