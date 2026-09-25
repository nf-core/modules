process VEPYR_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/da/dae579b8f2a3713d997875e16195f386b0fced234c8c9338549d16ed06eb9d83/data'
        : 'community.wave.seqera.io/library/htslib_vepyr:84d01ceaf76003ed'}"

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
    def reference = fasta ? "--fasta ${fasta}" : ''
    def version_arg = cache_version ? "--cache_version ${cache_version}" : ''
    def plugin_arg = plugin_cache ? "--plugin_cache_root ${plugin_cache}" : ''
    // --fork above 1 needs a tabix/CSI index: build one when none is given, and
    // run a single pipeline if the input cannot be indexed (not bgzip).
    def fork_command = tbi
        ? "fork=${task.cpus}"
        : "${vcf}".endsWith('.gz') ? "fork=${task.cpus}; tabix -f -p vcf ${vcf} || fork=1" : "fork=1"
    """
    ${fork_command}

    vepyr annotate \\
        -i ${vcf} \\
        -o ${prefix}.vcf.gz \\
        --dir_cache ${cache} \\
        ${reference} \\
        ${version_arg} \\
        ${plugin_arg} \\
        ${args} \\
        --fork \$fork \\
        --no_progress

    tabix ${args2} ${prefix}.vcf.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_vepyr"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
