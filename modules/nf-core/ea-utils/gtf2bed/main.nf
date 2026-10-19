process EAUTILS_GTF2BED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d2/d2a1f0d721596770e5289fe6369275008ffb9f3f23980bfb43fbf8d3d584a82e/data' :
        'community.wave.seqera.io/library/ea-utils:1.1.2.779--a489d71baa4a71a3' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("${prefix}.bed") , emit: bed
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    args   = task.ext.args ?: ''

    """
    echo $args
    """

    template 'gtf2bed.pl'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -e 'print substr(\$^V, 1)')
    END_VERSIONS
    """
}
