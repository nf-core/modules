process ELPREP_FILTER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/elprep:5.1.3--he881be0_1':
        'biocontainers/elprep:5.1.3--he881be0_1' }"

    input:
    tuple val(meta), path(bam), path(bai), path(target_regions_bed), path(filter_regions_bed), path(intermediate_bqsr_tables), path(recall_file)
    tuple val(meta2), path(reference_sequences)
    tuple val(meta3), path(reference_elfasta)
    tuple val(meta4), path(known_sites_elsites)
    val(run_haplotypecaller)
    val(run_bqsr)
    val(bqsr_tables_only)
    val(get_activity_profile)
    val(get_assembly_regions)


    output:
    tuple val(meta), path("*.{bam,sam}")            , emit: bam
    tuple val(meta), path("*.log")                  , emit: logs
    tuple val(meta), path("*.metrics.txt")          , optional: true, emit: metrics
    tuple val(meta), path("*.recall")               , optional: true, emit: recall
    tuple val(meta), path("*.vcf.gz")               , optional: true, emit: gvcf
    tuple val(meta), path("*.table")                , optional: true, emit: table
    tuple val(meta), path("*.activity_profile.igv") , optional: true, emit: activity_profile
    tuple val(meta), path("*.assembly_regions.igv") , optional: true, emit: assembly_regions
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains("--output-type sam") ? "sam" : "bam"

    // filter args
    def reference_sequences_cmd = reference_sequences ? "--replace-reference-sequences ${reference_sequences}" : ""
    def filter_regions_cmd      = filter_regions_bed  ? "--filter-non-overlapping-reads ${filter_regions_bed}" : ""

    // markdup args
    def markdup_cmd = args.contains("--mark-duplicates") ? "--mark-optical-duplicates ${prefix}.metrics.txt": ""

    // variant calling args
    def haplotyper_cmd = run_haplotypecaller ? "--haplotypecaller ${prefix}.g.vcf.gz": ""

    def fasta_cmd           = reference_elfasta   ? "--reference ${reference_elfasta}": ""
    def known_sites_cmd     = known_sites_elsites ? "--known-sites ${known_sites_elsites}": ""
    def target_regions_cmd  = target_regions_bed  ? "--target-regions ${target_regions_bed}": ""

    // bqsr args
    def bqsr_cmd = run_bqsr ? "--bqsr ${prefix}.recall": ""
    def bqsr_tables_only_cmd = bqsr_tables_only ? "--bqsr-tables-only ${prefix}.table": ""

    def intermediate_bqsr_cmd = intermediate_bqsr_tables ? "--bqsr-apply .": ""
    def input_recall_cmd = recall_file ? "--recal-file $recall_file" : ""
    // misc
    def activity_profile_cmd = get_activity_profile ? "--activity-profile ${prefix}.activity_profile.igv": ""
    def assembly_regions_cmd = get_assembly_regions ? "--assembly-regions ${prefix}.assembly_regions.igv": ""

    if ("$bam" == "${prefix}.${suffix}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    elprep filter ${bam} ${prefix}.${suffix} \\
        ${reference_sequences_cmd} \\
        ${filter_regions_cmd} \\
        ${markdup_cmd} \\
        ${haplotyper_cmd} \\
        ${fasta_cmd} \\
        ${known_sites_cmd} \\
        ${target_regions_cmd} \\
        ${bqsr_cmd} \\
        ${bqsr_tables_only_cmd} \\
        ${intermediate_bqsr_cmd} \\
        ${input_recall_cmd} \\
        ${activity_profile_cmd} \\
        ${assembly_regions_cmd} \\
        --nr-of-threads ${task.cpus} \\
        --log-path ./ \\
        $args

    mv logs/elprep/*.log .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        elprep: \$(elprep 2>&1 | head -n2 | tail -n1 |sed 's/^.*version //;s/ compiled.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains("--output-type sam") ? "sam" : "bam"
    def timestamp = "${java.time.OffsetDateTime.now().format(java.time.format.DateTimeFormatter.ISO_DATE_TIME)}"
    def markdup_cmd = args.contains("--mark-duplicates") ? "touch ${prefix}.metrics.txt": ""
    def bqsr_cmd = run_bqsr ? "touch ${prefix}.recall": ""
    def haplotyper_cmd = run_haplotypecaller ? "echo | gzip > ${prefix}.g.vcf.gz": ""
    def bqsr_tables_only_cmd = bqsr_tables_only ? "echo | gzip > ${prefix}.table": ""
    def activity_profile_cmd = get_activity_profile ? "touch ${prefix}.activity_profile.igv": ""
    def assembly_regions_cmd = get_assembly_regions ? "touch ${prefix}.assembly_regions.igv": ""

    if ("$bam" == "${prefix}.${suffix}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    touch ${prefix}.${suffix}
    touch elprep-${timestamp}.log
    ${markdup_cmd}
    ${bqsr_cmd}
    ${haplotyper_cmd}
    ${bqsr_tables_only_cmd}
    ${activity_profile_cmd}
    ${assembly_regions_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        elprep: \$(elprep 2>&1 | head -n2 | tail -n1 |sed 's/^.*version //;s/ compiled.*\$//')
    END_VERSIONS
    """
}
