process ELPREP_FILTER {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::elprep=5.1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/elprep:5.1.2--he881be0_0':
        'quay.io/biocontainers/elprep:5.1.2--he881be0_0' }"

    input:
    tuple val(meta), path(bam)
    val(run_haplotypecaller)
    val(run_bqsr)
    path(reference_sequences)
    path(filter_regions_bed)
    path(reference_elfasta)
    path(known_sites_elsites)
    path(target_regions_bed)
    path(intermediate_bqsr_tables)
    val(bqsr_tables_only)
    val(get_activity_profile)
    val(get_activity_regions)


    output:
    tuple val(meta), path("**.{bam,sam}")           ,emit: bam
    tuple val(meta), path("*.metrics.txt")          ,optional: true, emit: metrics
    tuple val(meta), path("*.recall")               ,optional: true, emit: recall
    tuple val(meta), path("*.vcf.gz")               ,optional: true, emit: gvcf
    tuple val(meta), path("*.table")                ,optional: true, emit: table
    tuple val(meta), path("*.activity_profile.igv") ,optional: true, emit: activity_profile
    tuple val(meta), path("*.activity_regions.igv") ,optional: true, emit: activity_regions
    path "versions.yml"                             ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains("--output-type sam") ? "sam" : "bam"

    // filter args
    if (reference_sequences) {
        args += " --replace-reference-sequences ${reference_sequences}"
    }
    if (filter_regions_bed) {
        args += " --filter-non-overlapping-reads ${filter_regions_bed}"
    }

    // markdup args
    if (args.contains("--mark-duplicates")){
        args += " --mark-optical-duplicates ${prefix}.metrics.txt"
    }
    // variant calling args
    if (run_haplotypecaller) {
        args += "--haplotypecaller ${prefix}.g.vcf.gz"
    }
    if (reference_elfasta) {
        args += " --reference ${reference_elfasta}"
    }
    if (known_sites_elsites) {
        args += " --known-sites ${known_sites_elsites}"
    }
    if (target_regions_bed) {
        args += " --target-regions ${target_regions_bed}"
    }
    // bqsr args
    if (run_bqsr) {
        args += " --bqsr ${prefix}.recall"
    }
    if (intermediate_bqsr_tables) {
        args += " --bqsr-apply ."
    }
    if (bqsr_tables_only) {
        args += " --bqsr-tables-only ${prefix}.table"
    }

    // misc
    if (get_activity_profile) {
        args += " ----activity-profile ${prefix}.activity_profile.igv"
    }
    if (get_activity_regions) {
        args += " --activity-regions ${prefix}.activity_regions.igv"
    }
    // errors
    if ((args.contains("--bqsr") || args.contains("--haplotypecaller")) && reference_elfasta == null) error "Reference required BQSR and variant calling."

    """
    elprep filter ${bam} ${prefix}.${suffix} \\
        --nr-of-threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        elprep: \$(elprep 2>&1 | head -n2 | tail -n1 |sed 's/^.*version //;s/ compiled.*\$//')
    END_VERSIONS
    """
}
