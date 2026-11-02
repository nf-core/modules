// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ANNOSINE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/71/71a65e004223b2b3d1acf1c4a0d8869b613c657c66c7f44454347715fde16202/data'
        : 'community.wave.seqera.io/library/annosine2:2.0.8--fef67f68ff78061f'}"

    input:
    tuple val(meta), path(fasta)
    val mode

    output:
    tuple val(meta), path("${prefix}.log"), emit: log
    tuple val(meta), path("${prefix}.fa") , emit: fa, optional: true
    tuple val("${task.process}"), val('annosine'), eval("python -c 'import importlib.metadata as m; print(m.version(\"annosine2\"))'"), emit: versions_annosine, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}_annosine"
    if ( "${fasta}" == "${prefix}.fa" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    AnnoSINE_v2 \\
        ${args} \\
        --threads ${task.cpus} \\
        ${mode} \\
        ${fasta} \\
        ${prefix} \\
        &> >(tee ${prefix}.log 2>&1)

    mv \\
        ${prefix}/Seed_SINE.fa \\
        ${prefix}.fa \\
        || echo 'AnnoSINE_v2 did not find SINE sequences. See log for details!'
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}_annosine"
    if ( "${fasta}" == "${prefix}.fa" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.log
    touch ${prefix}.fa
    """
}
