// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CUSTOM_SRATOOLSNCBISETTINGS {
    tag 'ncbi-settings'
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/12/1296a3f79e700d7479f4a294669da501d851daf93b0a246a4d951d630d908a53/data'
        : 'community.wave.seqera.io/library/sra-tools:3.2.1--c221ea63b92c74bf' }"

    input:
    val ids

    output:
    path('*.mkfg'), emit: ncbi_settings
    tuple val("${task.process}"), val('sratools'), eval("prefetch --version 2>&1 | grep -Eo '[0-9.]+'"), topic: versions, emit: versions_sratools

    when:
    task.ext.when == null || task.ext.when

    script:
    config = "/LIBS/GUID = \\\"${UUID.randomUUID().toString()}\\\"\\n/libs/cloud/report_instance_identity = \\\"true\\\"\\n"

    """
    echo ${config}
    """

    template "detect_ncbi_settings.sh"

    stub:
    """
    touch user-settings.mkfg
    """
}
