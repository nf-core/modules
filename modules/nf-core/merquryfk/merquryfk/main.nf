process MERQURYFK_MERQURYFK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fb/fb9dbdf098779f2b84401de5e6e8569dd9120018a47b92126876260a0bb621a6/data'
        : 'community.wave.seqera.io/library/merquryfk:1.2--5ac50fb470d146e0' }"

    input:
    tuple val(meta) , path(fastk_hist), path(fastk_ktab), path(assembly), path(haplotigs, arity:'0..*')
    tuple val(meta2), path(mathaptab) // optional, trio mode
    tuple val(meta3), path(pathaptab) // optional, trio mode

    output:
    tuple val(meta), path("${prefix}.completeness.stats")         , emit: stats
    tuple val(meta), path("${prefix}.*_only.bed")                 , emit: bed
    tuple val(meta), path("${prefix}.*.qv")                       , emit: assembly_qv
    tuple val(meta), path("${prefix}.qv")                         , emit: qv
    tuple val(meta), path("${prefix}.phased_block.bed")           , emit: phased_block_bed  , optional: true
    tuple val(meta), path("${prefix}.phased_block.stats")         , emit: phased_block_stats, optional: true
    tuple val(meta), path("*.{pdf,png}")                          , emit: images, optional: true
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('merquryfk'), val('1.2'), emit: versions_merquryfk, topic: versions
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('fastk'), val('1.2'), emit: versions_fastk, topic: versions
    tuple val("${task.process}"), val('R'), eval('R --version | sed "1!d; s/.*version //; s/ .*//"'), emit: versions_r, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if([mathaptab, pathaptab].any() && ![mathaptab, pathaptab].every()) {
        log.error("Error: Only one of the maternal and paternal hap tabs have been provided!")
    }

    def args        = task.ext.args ?: ''
    prefix          = task.ext.prefix ?: "${meta.id}"
    def fk_ktab     = fastk_ktab ? "${fastk_ktab.find { path -> path.toString().endsWith(".ktab") }}" : ''
    def mat_hapktab = mathaptab  ? "${mathaptab.find { path -> path.toString().endsWith(".ktab") }}"  : ''
    def pat_hapktab = pathaptab  ? "${pathaptab.find { path -> path.toString().endsWith(".ktab") }}"  : ''

    // Concatenate multiple haps if provided into a "pseudo" second haplotype for polyploid genome assessment
    def concatenate_haps = haplotigs.size() > 1
    if(concatenate_haps && haplotigs.collect { hap -> hap.getExtension() }.unique().size() > 1) {
        error("Error: Multiple haplotigs with different extensions found!")
    }
    def concat_out = haplotigs ? "temp.${haplotigs[0].getExtension() == 'gz' ? haplotigs[0].getName().split('\\.')[-2..-1].join('.') : haplotigs[0].getExtension()}" : ""
    def concat_command = concatenate_haps ? "cat ${haplotigs} > ${concat_out}" : ""
    def cleanup_command = concatenate_haps ? "rm ${concat_out}" : ""
    def hap_input = haplotigs && concatenate_haps ? "${concat_out}" : haplotigs
    """
    ${concat_command}

    MerquryFK \\
        ${args} \\
        -T$task.cpus \\
        ${fk_ktab} \\
        ${mat_hapktab} \\
        ${pat_hapktab} \\
        ${assembly} \\
        ${hap_input} \\
        ${prefix}

    ${cleanup_command}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.completeness.stats
    touch ${prefix}.qv
    touch ${prefix}._.qv
    touch ${prefix}._only.bed
    """
}
