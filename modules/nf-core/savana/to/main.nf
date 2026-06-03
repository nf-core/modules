process SAVANA_TO {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/savana:1.3.7--pyhdfd78af_0' :
        'quay.io/biocontainers/savana:1.3.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(tumour), path(tumour_index), path(snp_vcf), path(allele_counts_het_snps), path(breakpoints)
    tuple val(meta2), path(ref), path(ref_index)
    path(contigs)
    path(blacklist)
    val(g1000_vcf)


    output:
    tuple val(meta), path("${prefix}.sv_breakpoints.vcf"), emit: sv_breakpoints_vcf // savana run
    tuple val(meta), path("${prefix}.sv_breakpoints.bedpe"), emit: sv_breakpoints_bedpe // savana run
    tuple val(meta), path("${prefix}.sv_breakpoints_read_support.tsv"), emit: sv_breakpoints_read_support // savana run
    tuple val(meta), path("${prefix}.inserted_sequences.fa"), emit: inserted_sequences // savana run
    tuple val(meta), path("${prefix}.classified.vcf"), emit: classified_vcf // savana classify
    tuple val(meta), path("${prefix}.classified.somatic.vcf"), emit: somatic_vcf, optional: true // savana classify
    tuple val(meta), path("${prefix}.classified.somatic.bedpe"), emit: somatic_bedpe, optional: true // savana classify
    tuple val(meta), path("${prefix}.classified.germline.vcf"), emit: germline_vcf, optional: true // savana classify
    tuple val(meta), path("${prefix}.classified.{strict,lenient}.vcf"), emit: legacy_vcfs, optional: true // savana classify
    tuple val(meta), path("${prefix}.classified*.cna_rescue.vcf"), emit: cna_rescue_vcfs, optional: true
    tuple val(meta), path("${prefix}.somatic.labelled.vcf"), emit: labelled_vcf, optional: true // savana evaluate
    tuple val(meta), path("${prefix}.somatic.evaluation.stats"), emit: evaluation_stats, optional: true // savana evaluate
    tuple val(meta), path("${prefix}_allele_counts_hetSNPs.bed"), emit: allele_counts, optional: true
    tuple val(meta), path("${prefix}_raw_read_counts.tsv"), emit: raw_read_counts, optional: true
    tuple val(meta), path("${prefix}_read_counts_*_segmented.tsv"), emit: segmented_log2r, optional: true
    tuple val(meta), path("${prefix}_ranked_solutions.tsv"), emit: ranked_solutions, optional: true
    tuple val(meta), path("${prefix}_fitted_purity_ploidy.tsv"), emit: fitted_purity_ploidy, optional: true
    tuple val(meta), path("${prefix}_segmented_absolute_copy_number.tsv"), emit: cna, optional: true
    tuple val(meta), path("*kbp_bin_ref_*_${prefix}*.bed"), emit: binned_ref, optional: true
    tuple val("${task.process}"), val('savana'), eval("python -c \"import importlib.metadata as m; print(m.version('savana'))\""), emit: versions_savana, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def contigs_arg = contigs ? "--contigs ${contigs}" : ""
    def allele_arg = snp_vcf
        ? "--snp_vcf ${snp_vcf}"
        : g1000_vcf
            ? "--g1000_vcf ${g1000_vcf}"
            : allele_counts_het_snps ? "--allele_counts_het_snps ${allele_counts_het_snps}" : ""
    def blacklist_arg = blacklist ? "--blacklist ${blacklist}" : ""
    def ref_index_arg = ref_index ? "--ref_index ${ref_index}" : ""
    def breakpoints_arg = breakpoints ? "--breakpoints ${breakpoints}" : ""

    """
    savana to \\
        --tumour ${tumour} \\
        --ref ${ref} \\
        ${ref_index_arg} \\
        --outdir "./outdir" \\
        --sample ${prefix} \\
        --threads ${task.cpus} \\
        --cna_threads ${task.cpus} \\
        ${contigs_arg} \\
        ${blacklist_arg} \\
        ${allele_arg} \\
        ${breakpoints_arg} \\
        ${args}

    mv ./outdir/* .
    rmdir ./outdir
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def cna_outputs = (snp_vcf || g1000_vcf || allele_counts_het_snps) ? """
    touch ${prefix}_raw_read_counts.tsv
    touch ${prefix}_read_counts_self_log2r_segmented.tsv
    touch ${prefix}_ranked_solutions.tsv
    touch ${prefix}_fitted_purity_ploidy.tsv
    touch ${prefix}_segmented_absolute_copy_number.tsv
    touch ${prefix}_allele_counts_hetSNPs.bed
    touch 10kbp_bin_ref_all_${prefix}_with_SV_breakpoints.bed
    """ : ""

    """
    touch ${prefix}.sv_breakpoints.vcf
    touch ${prefix}.sv_breakpoints.bedpe
    touch ${prefix}.sv_breakpoints_read_support.tsv
    touch ${prefix}.inserted_sequences.fa
    touch ${prefix}.classified.vcf
    touch ${prefix}.classified.somatic.vcf
    touch ${prefix}.classified.somatic.bedpe
    touch ${prefix}.classified.strict.vcf
    touch ${prefix}.classified.lenient.vcf
    ${cna_outputs}
    """
}
