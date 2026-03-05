process BAMCMP {
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6e/6e5ee3676abe7e65f65eca55e8dbc76f4dd195a44679cd1a785943d7a0d598f1/data' :
        'community.wave.seqera.io/library/bamcmp_samtools:2f211ea999bb54f5' }"

    input:
    tuple val(meta), path(primary_aligned_bam), path(contaminant_aligned_bam)

    output:
    tuple val(meta), path("${prefix}.bam") , emit: primary_filtered_bam
    tuple val(meta), path("${prefix2}.bam"), emit: contamination_bam
    tuple val("${task.process}"), val('bamcmp'), eval('echo 2.2'), topic: versions, emit: versions_bamcmp
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_primary"
    prefix2 = task.ext.prefix2 ?: "${meta.id}_contaminant"

    if ("$primary_aligned_bam" == "${prefix}.bam"  | "$contaminant_aligned_bam" == "${prefix}.bam"  ) {
        error "Input and output names for the primary-genome bam file are the same, use \"task.ext.prefix\" to disambiguate!"
    }
    if ("$primary_aligned_bam" == "${prefix2}.bam" | "$contaminant_aligned_bam" == "${prefix2}.bam" ) {
        error "Input and output names for the contaminant-genome bam file are the same, use \"task.ext.prefix2\" to disambiguate!"
    }
    if ("primary_aligned_bam"    == "contaminant_aligned_bam" ) {
        error "Input file names for the two bam files are the same, ensure they are renamed upstream."
    }
    if ("${prefix}.bam"    == "${prefix2}.bam" ) {
        error "Output names for the two bam files are identical, use \"task.ext.prefix\" and \"task.ext.prefix2\" to disambiguate!"
    }
    """
    bamcmp \\
        -1 $primary_aligned_bam \\
        -2 $contaminant_aligned_bam \\
        -A ${prefix}.sam \\
        -B ${prefix2}.sam \\
        $args

    # Convert SAM outputs to BAM to ensure proper BAM format
    samtools view -b -o ${prefix}.bam ${prefix}.sam
    samtools view -b -o ${prefix2}.bam ${prefix2}.sam
    rm ${prefix}.sam ${prefix2}.sam
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_primary"
    prefix2 = task.ext.prefix2 ?: "${meta.id}_contaminant"

    if ("$primary_aligned_bam" == "${prefix}.bam"  | "$contaminant_aligned_bam" == "${prefix}.bam"  )
        error "Input and output names for the primary-genome bam file are the same, use \"task.ext.prefix\" to disambiguate!"
    if ("$primary_aligned_bam" == "${prefix2}.bam" | "$contaminant_aligned_bam" == "${prefix2}.bam" )
        error "Input and output names for the contaminant-genome bam file are the same, use \"task.ext.prefix2\" to disambiguate!"
    if ("primary_aligned_bam"    == "contaminant_aligned_bam" )
        error "Input file names for the two bam files are the same, ensure they are renamed upstream."
    if ("${prefix}.bam"    == "${prefix2}.bam" )
        error "Output names for the two bam files are identical, use \"task.ext.prefix\" and \"task.ext.prefix2\" to disambiguate!"
    """
    touch ${prefix}.bam
    touch ${prefix2}.bam
    """

}
