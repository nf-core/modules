// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process METAPHLAN_MAKEDB {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1a/1a05e0cb9a7387c4711a0bed4e283152bba02d06712d97c0253aa481b48e1a3e/data'
        : 'community.wave.seqera.io/library/metaphlan:4.3--6cd977ac2f320230'}"

    input:
    val meta

    output:
    tuple val(meta), path("metaphlan_db_latest"), emit: db
    tuple val("${task.process}"), val('metaphlan'), eval("metaphlan --version 2>&1 | head -n 1 | cut -d ' ' -f 3"), emit: versions_metaphlan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    metaphlan \\
        --install \\
        --nproc $task.cpus \\
        --db_dir metaphlan_db_latest \\
        ${args}

    """

    stub:
    def args = task.ext.args ?: ''
    def mmi = args.contains('--long_reads') ? 'touch metaphlan_db_latest/mpa_vJan26_CHOCOPhlAnSGB_202611_xasm20_B3_O3_12.mmi' : ''
    """
    echo "${args}"
    mkdir metaphlan_db_latest
    touch metaphlan_db_latest/mpa_latest
    touch metaphlan_db_latest/mpa_vJan26_CHOCOPhlAnSGB_202611.1.bt2l
    touch metaphlan_db_latest/mpa_vJan26_CHOCOPhlAnSGB_202611.2.bt2l
    touch metaphlan_db_latest/mpa_vJan26_CHOCOPhlAnSGB_202611.3.bt2l
    touch metaphlan_db_latest/mpa_vJan26_CHOCOPhlAnSGB_202611.4.bt2l
    touch metaphlan_db_latest/mpa_vJan26_CHOCOPhlAnSGB_202611.rev.1.bt2l
    touch metaphlan_db_latest/mpa_vJan26_CHOCOPhlAnSGB_202611.rev.2.bt2l
    touch metaphlan_db_latest/mpa_vJan26_CHOCOPhlAnSGB_202611.fna
    touch metaphlan_db_latest/mpa_vJan26_CHOCOPhlAnSGB_202611.nwk
    touch metaphlan_db_latest/mpa_vJan26_CHOCOPhlAnSGB_202611.pkl
    touch metaphlan_db_latest/mpa_vJan26_CHOCOPhlAnSGB_202611_VINFO.csv
    touch metaphlan_db_latest/mpa_vJan26_CHOCOPhlAnSGB_202611_VSG.fna
    ${mmi}

    """
}
