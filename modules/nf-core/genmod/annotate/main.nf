// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process GENMOD_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ff/ffda9280a74a528145dfb7125bfaab17be55aefd0860d56225c37e2a9d0964ec/data' :
        'community.wave.seqera.io/library/genmod:3.13.0--abcc76f3e3411452' }"

    input:
    tuple val(meta), path(input_vcf)

    output:
    tuple val(meta), path("*_annotate.vcf"), emit: vcf
    tuple val("${task.process}"), val('genmod'), eval("genmod --version | sed 's/^.*genmod version: //'"), emit: versions_genmod, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    genmod \\
        annotate \\
        $args \\
        --outfile ${prefix}_annotate.vcf \\
        $input_vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_annotate.vcf
    """
}
