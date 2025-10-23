process EIDO_CONVERT {
    tag "$samplesheet"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/eido_peppy:2522b1352d5d6547' :
        'community.wave.seqera.io/library/eido_peppy:3721c3f85cc3d076' }"

    input:
    path samplesheet
    val format

    output:
    path "versions.yml"           , emit: versions
    path "${prefix}.${format}"    , emit: samplesheet_converted

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "samplesheet_converted"
    """
    eido \\
        convert \\
        -f $format \\
        $samplesheet \\
        $args \\
        -p samples=${prefix}.${format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eido: \$(echo \$(eido --version 2>&1) | sed 's/^.*eido //;s/ .*//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "samplesheet_converted"
    """
    touch ${prefix}.${format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eido: \$(echo \$(eido --version 2>&1) | sed 's/^.*eido //;s/ .*//' )
    END_VERSIONS
    """
}
