// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process STRDROP_CALL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/99/9915c87a51b3c5af344c3044ac79c7badb8b44a1aed0eeba4363bc10e2a4bc95/data'
        : 'community.wave.seqera.io/library/pip_python_typing_extensions_numpy_strdrop:de2d4e2fcc806913'}"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(training_set_json)
    tuple val(meta3), path(training_set_vcfs, stageAs: 'input/*')

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('strdrop'), eval("strdrop --version | sed 's/.* //g'"), topic: versions, emit: versions_strdrop

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def training_set = training_set_json ? "--training-set ${training_set_json}" : '--training-set ./input'

    if (training_set_json && training_set_vcfs) {
        error("Please provide only one of 'training_set_json' or 'training_set_vcfs' as training set input.")
    }
    """
    strdrop \\
        call \\
        $args \\
        $training_set \\
        $vcf \\
        ${prefix}.vcf.gz
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (training_set_json && training_set_vcfs) {
        error("Please provide only one of 'training_set_json' or 'training_set_vcfs' as training set input.")
    }
    """
    echo $args

    echo "" | gzip > ${prefix}.vcf.gz
    """
}
