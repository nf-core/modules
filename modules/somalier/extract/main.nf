
process SOMALIER_EXTRACT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "YOUR-TOOL-HERE" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
          'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
          'brentp/somalier:latest' }"


    input:
    tuple val(meta), path(bam), path(bai)
    path(ref)
    path(refidx)
    path(sites)

	output:
    tuple val(meta), path("*.somalier.extract/*.somalier"),             emit: extract
    path "versions.yml",                                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    somalier extract ${args} -d ${prefix}.somalier.extract \
    --sites ${sites} \
    -f ${ref} \
    ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(echo \$(somalier_dbg4 2>&1) | sed 's/^.*somalier version: //; s/Commands:.*\$//')
    END_VERSIONS
    """

}
