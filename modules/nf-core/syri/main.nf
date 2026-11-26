// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SYRI {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/44/444b93baa5cb545d52452da19ff2516424d2f14513ad9698b9933152aad06b39/data'
        : 'community.wave.seqera.io/library/syri_pandas:84f5dfbb480b2430'}"

    input:
    tuple val(meta), path(infile)
    tuple val(meta2), path(query_fasta)
    tuple val(meta3), path(reference_fasta)
    val(file_type)

    output:
    tuple val(meta), path("*.syri.out"), emit: syri , optional: true
    tuple val(meta), path("*.error.log"), emit: error , optional: true
    tuple val("${task.process}"), val('syri'), eval("syri --version"), emit: versions_syri, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ( ! ( file_type in [ 'T', 'S', 'B', 'P' ] ) ) { error "File type should be one of [ 'T', 'S', 'B', 'P' ]" }
    """
    syri \\
        -c ${infile} \\
        -q ${query_fasta} \\
        -r ${reference_fasta} \\
        -F ${file_type} \\
        ${args} \\
        --prefix ${prefix} \\
        2>| >(tee "${prefix}.error.log" >&2) \\
        || echo "Errors from syri printed to ${prefix}.error.log"

    [ -f "${prefix}.syri.out" ] \\
        && rm "${prefix}.error.log" \\
        || echo 'Syri failed and no syri.out file was created'
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.syri.out
    """
}
