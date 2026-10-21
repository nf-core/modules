process CUSTOM_GETCHROMSIZES {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fc/fca45dc02d0b0906faa49dcbb706dc4dc07b2327b1884e8b344893d10fc87ae1/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.22.1--ef1e734834ad9fa5'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path ("*.sizes"), emit: sizes
    tuple val(meta), path ("*.fai")  , emit: fai
    tuple val(meta), path ("*.gzi")  , emit: gzi, optional: true
    path  "versions.yml"             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
WARNING: The getchromsizes process has been deprecated. Please use nf-core/modules/samtools/faidx.

Reason:
Getting chromosome sizes was added to samtools/faidx (https://github.com/nf-core/modules/pull/7041)
via a boolean switch, making 'getchromsizes' unnecessary.
"""
    assert false: deprecation_message

    """
    samtools faidx $fasta
    cut -f 1,2 ${fasta}.fai > ${fasta}.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getchromsizes: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def deprecation_message = """
WARNING: The getchromsizes process has been deprecated. Please use nf-core/modules/samtools/faidx.

Reason:
Getting chromosome sizes was added to samtools/faidx (https://github.com/nf-core/modules/pull/7041)
via a boolean switch, making 'getchromsizes' unnecessary.
"""
    assert false: deprecation_message
}
