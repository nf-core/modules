process VEPYR_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/da/dae579b8f2a3713d997875e16195f386b0fced234c8c9338549d16ed06eb9d83/data'
        : 'community.wave.seqera.io/library/htslib_vepyr:84d01ceaf76003ed'}"

    input:
    // Staged under input/ so an input named ${prefix}.vcf.gz (e.g. the output of
    // an upstream module using the same meta.id) never collides with -o.
    tuple val(meta), path(vcf, stageAs: 'input/*'), path(tbi, stageAs: 'input/*')
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
    prefix = task.ext.prefix ?: "${meta.id}"
    // vepyr opens the reference through its .fai, and a bgzip FASTA through its
    // .gzi as well, and builds neither. fai and gzi are never named on the
    // command line: they are inputs so that they are staged alongside the FASTA,
    // without which --everything/--hgvsc fail.
    def reference = fasta ? "--fasta ${fasta}" : ''
    def version_arg = cache_version ? "--cache_version ${cache_version}" : ''
    def plugin_arg = plugin_cache ? "--plugin_cache_root ${plugin_cache}" : ''
    // --fork above 1 requires a tabix/CSI index on the input; vepyr raises
    // without one. Fall back to a single pipeline so a missing index costs
    // throughput rather than failing the task.
    //
    // The flag is emitted *after* ${args} rather than before it, the one place
    // this module overrides the user: --fork and --workers are a single
    // argparse option, so the last occurrence wins, and an ext.args carrying
    // either spelling would silently defeat the fallback and fail the task.
    // Parallelism belongs to the cpus directive here.
    def fork = tbi ? task.cpus : 1
    """
    vepyr annotate \\
        -i ${vcf} \\
        -o ${prefix}.vcf.gz \\
        --dir_cache ${cache} \\
        ${reference} \\
        ${version_arg} \\
        ${plugin_arg} \\
        ${args} \\
        --fork ${fork} \\
        --no_progress

    tabix ${args2} ${prefix}.vcf.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
