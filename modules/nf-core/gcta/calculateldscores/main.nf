process GCTA_CALCULATELDSCORES {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed2d6a4b4f43f3230014aa67bb40feb0acbd7cc41ef0b0a895002f1befc6502c/data'
        : 'community.wave.seqera.io/library/gcta_r-base:31127c93877b38de'}"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    val ld_score_region

    output:
    tuple val(meta), path("*_gcta_ld.score.ld"), emit: ld_scores
    tuple val(meta), path("*_snp_group*.txt"), emit: snp_group_files
    tuple val("${task.process}"), val("gcta"), eval("gcta --version | sed -En 's/^[*] version v([0-9.]*).*/\\1/p'"), emit: versions_gcta, topic: versions
    tuple val("${task.process}"), val("r-base"), eval('Rscript -e "cat(as.character(getRversion()))"'), emit: versions_rbase, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bfile_prefix = bed.baseName
    """
    gcta \\
        --bfile ${bfile_prefix} \\
        --ld-score-region ${ld_score_region} \\
        --out ${prefix}_gcta_ld \\
        --thread-num ${task.cpus} \\
        ${args}

    Rscript - <<'EOF'
    lds_seg <- read.table(
      "${prefix}_gcta_ld.score.ld",
      header = TRUE,
      colClasses = c("character", rep("numeric", 8))
    )

    ld_score_quantiles <- quantile(lds_seg\$ldscore_SNP, probs = c(0.25, 0.50, 0.75), names = FALSE)

    snp_group_filters <- list(
      lds_seg\$ldscore_SNP <= ld_score_quantiles[1],
      lds_seg\$ldscore_SNP > ld_score_quantiles[1] & lds_seg\$ldscore_SNP <= ld_score_quantiles[2],
      lds_seg\$ldscore_SNP > ld_score_quantiles[2] & lds_seg\$ldscore_SNP <= ld_score_quantiles[3],
      lds_seg\$ldscore_SNP > ld_score_quantiles[3]
    )

    for (group_index in seq_along(snp_group_filters)) {
      writeLines(
        lds_seg\$SNP[snp_group_filters[[group_index]]],
        sprintf("${prefix}_snp_group%d.txt", group_index)
      )
    }
    EOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "%s\\n" \
        "SNP A1 A2 Freq b se p N ldscore_SNP" \
        "stub_snp1 A G 0.10 0.01 0.02 0.50 100 1.00" \
        "stub_snp2 A G 0.20 0.02 0.03 0.40 100 2.00" \
        "stub_snp3 A G 0.30 0.03 0.04 0.30 100 3.00" \
        "stub_snp4 A G 0.40 0.04 0.05 0.20 100 4.00" \
        > ${prefix}_gcta_ld.score.ld
    printf "stub_snp1\n" > ${prefix}_snp_group1.txt
    printf "stub_snp2\n" > ${prefix}_snp_group2.txt
    printf "stub_snp3\n" > ${prefix}_snp_group3.txt
    printf "stub_snp4\n" > ${prefix}_snp_group4.txt
    """
}
