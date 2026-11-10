// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ATLAS_CALL {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a4/a494a4d166b29fdf9a27ac525c7aa4dc8ec476d0e624d174e6c72a0e87499dab/data'
        : 'community.wave.seqera.io/library/atlas:0.9.9--18077f78d3bd1a61'}"

    input:
    tuple val(meta), path(bam), path(bai), path(pmd), path(recal)
    path fasta
    path fai
    path known_alleles
    val method

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('atlas'), eval("(atlas 2>&1) | grep Atlas | head -n 1 | sed -e 's/^[[:space:]]*Atlas //'"), emit: versions_atlas, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args   ?: ''
    def recal_file         = recal           ? "recal=${recal}" : ""
    def pmd_file           = pmd             ? "pmdFile=${pmd}" : ""
    def known_alleles_file = known_alleles   ? "alleles=${known_alleles}" : ""

    def valid_method = ['MLE', 'Bayesian', 'allelePresence', 'randomBase', 'majorityBase']
    if ( !valid_method.contains(method) )  { error "Unrecognised calling method for ATLAS_CALL. Options: MLE, Bayesian, allelePresence, randomBase, majorityBase" }

    """
    atlas \\
        task=call \\
        bam=${bam} \\
        fasta=${fasta} \\
        ${recal_file} \\
        ${pmd_file} \\
        ${known_alleles_file} \\
        method=${method} \\
        ${args}
    """

}
