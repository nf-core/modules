process MSISENSOR_SCAN {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/msisensor:0.5--hb3646a4_2'
        : 'biocontainers/msisensor:0.5--hb3646a4_2'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tab"), emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
        WARNING: This module has been deprecated. Please use nf-core/modules/nf-core/msisensor2/scan or nf-core/modules/nf-core/msisensorpro/scan

        Reason:
        This module is no longer fit for purpose as it is now archived and no longer maintained.
    """
    assert false: deprecation_message

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    msisensor \\
        scan \\
        -d ${fasta} \\
        -o ${prefix}.msisensor_scan.tab \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor: \$(msisensor 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p')
    END_VERSIONS
    """

    stub:
    def deprecation_message = """
        WARNING: This module has been deprecated. Please use nf-core/modules/nf-core/msisensor2/scan or nf-core/modules/nf-core/msisensorpro/scan

        Reason:
        This module is no longer fit for purpose as it is now archived and no longer maintained.
    """
    assert false: deprecation_message
}
