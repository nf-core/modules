process ENSEMBLVEP_FILTERVEP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4b/4b5a8c173dc9beaa93effec76b99687fc926b1bd7be47df5d6ce19d7d6b4d6b7/data'
        : 'community.wave.seqera.io/library/ensembl-vep:115.2--90ec797ecb088e9a'}"

    input:
    tuple val(meta), path(input)
    path feature_file

    output:
    tuple val(meta), path("*.${extension}"), emit: output
    tuple val("${task.process}"), val('ensemblvep'), eval("vep --help | sed -n '/ensembl-vep/s/.*: //p'"), topic: versions, emit: versions_ensemblvep

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    extension = task.ext.suffix ?: "vcf"
    """
    touch ${prefix}.${extension}
    """
}
