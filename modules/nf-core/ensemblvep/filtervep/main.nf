process ENSEMBLVEP_FILTERVEP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/84/8480a9566a1d6162c5119c7ec7e2a9991c6ed4415eabe3b11c728bed2a85b77f/data'
        : 'community.wave.seqera.io/library/ensembl-vep:115.1--dae34d4b1c5d54d2'}"

    input:
    tuple val(meta), path(input)
    path feature_file

    output:
    tuple val(meta), path("*.${extension}"), emit: output
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    extension = task.ext.suffix ?: "vcf"
    """
    filter_vep \\
        ${args} \\
        --input_file ${input} \\
        --output_file ${prefix}.${extension} \\
        --only_matched

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    extension = task.ext.suffix ?: "vcf"
    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}
