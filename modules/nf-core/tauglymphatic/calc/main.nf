// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process TAUGLYMPHATIC_CALC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/db/db230f0ed6fef58e6960f9a04a15307470749a848c165ea6c2964803d5505471/data'
        : 'community.wave.seqera.io/library/tau-glymphatic-continuum:1.1.2--19d94451d3401cc3'}"

    input:
    tuple val(meta), val(alps), val(apoe_e4), val(trem2_r47h), val(nlrp3)

    output:
    tuple val(meta), path("${prefix}.{txt,json}"), emit: report
    tuple val("${task.process}"), val('tau-glymphatic-continuum'), eval("python3 -c \"import importlib.metadata; print(importlib.metadata.version('tau-glymphatic-continuum'))\""), topic: versions, emit: versions_tauglymphatic

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def alps_arg = alps != null ? "--alps ${alps}" : ''
    def apoe_arg = apoe_e4 != null ? "--apoe-e4 ${apoe_e4}" : ''
    def trem2_arg = trem2_r47h ? '--trem2-r47h' : ''
    def nlrp3_arg = nlrp3 ? '--nlrp3' : ''
    def ext = (args.contains('--json') && !args.contains('--audit')) ? 'json' : 'txt'
    """
    tau-calc \\
        ${alps_arg} \\
        ${apoe_arg} \\
        ${trem2_arg} \\
        ${nlrp3_arg} \\
        ${args} \\
        > ${prefix}.${ext}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ext = (args.contains('--json') && !args.contains('--audit')) ? 'json' : 'txt'
    """
    touch ${prefix}.${ext}
    """
}
