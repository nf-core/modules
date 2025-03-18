process MERQURYFK_MERQURYFK {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/11/116bae9e8314713209c08c429a95bbf8ce77486872bede198bd86f8267cffb7b/data' :
        'community.wave.seqera.io/library/fastk_merquryfk_r-argparse_r-base_pruned:ec47d7677abb1c46' }"

    input:
    tuple val(meta) , path(fastk_hist),path(fastk_ktab), path(assembly), path(haplotigs)
    tuple val(meta2), path(mathaptab) // optional, trio mode
    tuple val(meta3), path(pathaptab) // optional, trio mode

    output:
    tuple val(meta), path("${prefix}.completeness.stats")         , emit: stats
    tuple val(meta), path("${prefix}.*_only.bed")                 , emit: bed
    tuple val(meta), path("${prefix}.*.qv")                       , emit: assembly_qv
    tuple val(meta), path("${prefix}.*.spectra-cn.fl.{png,pdf}")  , emit: spectra_cn_fl     , optional: true
    tuple val(meta), path("${prefix}.*.spectra-cn.ln.{png,pdf}")  , emit: spectra_cn_ln     , optional: true
    tuple val(meta), path("${prefix}.*.spectra-cn.st.{png,pdf}")  , emit: spectra_cn_st     , optional: true
    tuple val(meta), path("${prefix}.qv")                         , emit: qv
    tuple val(meta), path("${prefix}.spectra-asm.fl.{png,pdf}")   , emit: spectra_asm_fl    , optional: true
    tuple val(meta), path("${prefix}.spectra-asm.ln.{png,pdf}")   , emit: spectra_asm_ln    , optional: true
    tuple val(meta), path("${prefix}.spectra-asm.st.{png,pdf}")   , emit: spectra_asm_st    , optional: true
    tuple val(meta), path("${prefix}.phased_block.bed")           , emit: phased_block_bed  , optional: true
    tuple val(meta), path("${prefix}.phased_block.stats")         , emit: phased_block_stats, optional: true
    tuple val(meta), path("${prefix}.continuity.N.{pdf,png}")     , emit: continuity_N      , optional: true
    tuple val(meta), path("${prefix}.block.N.{pdf,png}")          , emit: block_N           , optional: true
    tuple val(meta), path("${prefix}.block.blob.{pdf,png}")       , emit: block_blob        , optional: true
    tuple val(meta), path("${prefix}.hapmers.blob.{pdf,png}")     , emit: hapmers_blob      , optional: true
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    if([mathaptab, pathaptab].any() && ![mathaptab, pathaptab].every()) {
        log.error("Error: Only one of the maternal and paternal hap tabs have been provided!")
    }

    fk_ktab     = fastk_ktab ? "${fastk_ktab.find{ it.toString().endsWith(".ktab") }}" : ''
    mat_hapktab = mathaptab  ? "${mathaptab.find{ it.toString().endsWith(".ktab") }}" : ''
    pat_hapktab = pathaptab  ? "${pathaptab.find{ it.toString().endsWith(".ktab") }}" : ''
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def FASTK_VERSION   = '1.1.0'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def MERQURY_VERSION = '1.1.1'
    """
    MerquryFK \\
        $args \\
        -T$task.cpus \\
        ${fk_ktab} \\
        ${mat_hapktab} \\
        ${pat_hapktab} \\
        $assembly \\
        $haplotigs \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
        merquryfk: $MERQURY_VERSION
        r: \$( R --version | sed '1!d; s/.*version //; s/ .*//' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def FASTK_VERSION   = '1.1.0'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def MERQURY_VERSION = '1.1.1'
    """
    touch ${prefix}.completeness.stats
    touch ${prefix}.qv
    touch ${prefix}._.qv
    touch ${prefix}._only.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
        merquryfk: $MERQURY_VERSION
        r: \$( R --version | sed '1!d; s/.*version //; s/ .*//' )
    END_VERSIONS
    """
}
