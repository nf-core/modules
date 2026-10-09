// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VEPYR_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bb/bbf0f7f5a0c6b6586735b7636b219e3b788b9c9c1a86994b776556c63892fb4b/data'
        : 'community.wave.seqera.io/library/htslib_vepyr:408f2021357958aa'}"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(cache)
    tuple val(meta3), path(fasta), path(fai), path(gzi)
    val cache_version
    tuple val(meta4), path(plugin_cache)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('vepyr'), eval("vepyr --version | cut -d' ' -f2"), topic: versions, emit: versions_vepyr
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'"), topic: versions, emit: versions_tabix

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_vepyr"
    if ("${vcf}" == "${prefix}.vcf.gz") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    if (!fasta) {
        error("VEPYR_ANNOTATE requires a reference FASTA: vepyr always annotates with --everything, which needs it.")
    }
    def version_arg = cache_version ? "--cache_version ${cache_version}" : ''
    def plugin_arg = plugin_cache ? "--plugin_cache_root ${plugin_cache}" : ''
    // --fork above 1 needs a tabix/CSI index, so a plain .vcf runs one pipeline
    def fork = tbi || "${vcf}".endsWith('.gz') ? task.cpus : 1
    def index_command = !tbi && "${vcf}".endsWith('.gz') ? "tabix -p vcf ${vcf}" : ''
    """
    ${index_command}

    vepyr annotate \\
        -i ${vcf} \\
        -o ${prefix}.vcf.gz \\
        --dir_cache ${cache} \\
        --fasta ${fasta} \\
        ${version_arg} \\
        ${plugin_arg} \\
        ${args} \\
        --fork ${fork} \\
        --no_progress

    tabix ${args2} ${prefix}.vcf.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_vepyr"
    if ("${vcf}" == "${prefix}.vcf.gz") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
