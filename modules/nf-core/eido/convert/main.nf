process EIDO_CONVERT {
    tag "$samplesheet"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4a/4aef51a3d75d6b950bb643ed5ebc1d7243d67bbf3de0410fcaa7d347e8fc0007/data' :
        'community.wave.seqera.io/library/eido_peppy_setuptools:53ed68799568c4fa' }"

    input:
    path samplesheet
    val format
    path pep_input_base_dir

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
