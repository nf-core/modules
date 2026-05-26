process SAVANA_CNA {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/savana:1.3.7--pyhdfd78af_0'
        : 'quay.io/biocontainers/savana:1.3.7--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(tumour), path(tumour_index), path(normal), path(normal_index), path(snp_vcf), path(allele_counts_het_snps), path(breakpoints)
    tuple val(meta2), path(ref), path(ref_index)
    path blacklist
    val g1000_vcf

    output:
    tuple val(meta), path("cna_out/*segmented_absolute_copy_number.tsv"), emit: cna
    tuple val(meta), path("cna_out/*ranked_solutions.tsv"), emit: ranked_solutions
    tuple val(meta), path("cna_out/*fitted_purity_ploidy.tsv"), emit: fitted_purity_ploidy
    tuple val(meta), path("cna_out/*read_counts_*_log2r_segmented.tsv"), emit: segmented_log2r
    tuple val(meta), path("cna_out/*raw_read_counts.tsv"), emit: raw_read_counts
    tuple val(meta), path("cna_out/*kbp_bin_ref_*_*.bed"), emit: binned_ref
    tuple val(meta), path("cna_out/*allele_counts_hetSNPs.bed"), optional: true, emit: allele_counts
    tuple val("${task.process}"), val("savana"), eval("python -c \"import importlib.metadata as m; print(m.version('savana'))\""), emit: versions_savana, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ""
    def normal_arg = normal ? "--normal ${normal}" : ""
    def allele_arg = snp_vcf
        ? "--snp_vcf ${snp_vcf}"
        : g1000_vcf
            ? "--g1000_vcf ${g1000_vcf}"
            : allele_counts_het_snps ? "--allele_counts_het_snps ${allele_counts_het_snps}" : ""
    def blacklist_arg = blacklist ? "--blacklist ${blacklist}" : ""
    def breakpoints_arg = breakpoints ? "--breakpoints ${breakpoints}" : ""
    """
    savana cna \
      --tumour ${tumour} \
      ${normal_arg} \
      --ref ${ref} \
      --outdir cna_out \
      --sample ${prefix} \
      --cna_threads ${task.cpus} \
      ${allele_arg} \
      ${blacklist_arg} \
      ${breakpoints_arg} \
      ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    def norm_mode = normal ? "mnorm" : "self"
    def do_allele = (snp_vcf || g1000_vcf || allele_counts_het_snps)

    """
    mkdir -p cna_out

    printf "bin\tchromosome\tstart\tend\tperc_known_bases\tuse_bin\ttumour_read_count\tnormal_read_count\n" > cna_out/${prefix}_raw_read_counts.tsv
    printf "bin\tchromosome\tstart\tend\tperc_known_bases\tuse_bin\tlog2r_copynumber\tseg_id\tseg_log2r_copynumber\n" > cna_out/${prefix}_read_counts_${norm_mode}_log2r_segmented.tsv

    touch cna_out/10kbp_bin_ref_all_${prefix}.bed

    printf "purity\tploidy\tdistance\trank\n" > cna_out/${prefix}_ranked_solutions.tsv
    printf "purity\tploidy\tdistance\trank\n" > cna_out/${prefix}_fitted_purity_ploidy.tsv

    if [ "${do_allele}" = "true" ]; then
      printf "chromosome\tstart\tend\tsegment_id\tbin_count\tsum_of_bin_lengths\tweight\tcopyNumber\tminorAlleleCopyNumber\tmeanBAF\tno_hetSNPs\n" > cna_out/${prefix}_segmented_absolute_copy_number.tsv
      touch cna_out/${prefix}_allele_counts_hetSNPs.bed
    else
      printf "chromosome\tstart\tend\tsegment_id\tbin_count\tsum_of_bin_lengths\tweight\tcopyNumber\n" > cna_out/${prefix}_segmented_absolute_copy_number.tsv
    fi
    """
}
