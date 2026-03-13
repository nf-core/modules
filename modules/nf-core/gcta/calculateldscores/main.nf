process GCTA_CALCULATELDSCORES {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/gcta_r-base:31127c93877b38de' :
        'community.wave.seqera.io/library/gcta_r-base:31127c93877b38de' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    val ld_score_region

    output:
    tuple val(meta), path("${meta.id}_gcta_ld.score.ld"), emit: ld_scores
    tuple val(meta), path("${meta.id}_snp_group1.txt"), path("${meta.id}_snp_group2.txt"), path("${meta.id}_snp_group3.txt"), path("${meta.id}_snp_group4.txt"), emit: snp_group_files
    tuple val("${task.process}"), val("gcta"), eval("gcta --version 2>&1 | head -n 1"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def extra_args = task.ext.args ?: ''

    """
    set -euo pipefail

    gcta \\
        --bfile ${meta.id} \\
        --ld-score-region ${ld_score_region} \\
        --out ${meta.id}_gcta_ld \\
        --thread-num ${task.cpus} ${extra_args}

    Rscript - ${meta.id}_gcta_ld.score.ld ${meta.id} <<'EOF'
    args <- commandArgs(trailingOnly = TRUE)
    filename <- args[1]
    out_prefix <- args[2]

    lds_seg <- read.table(
      filename,
      header = TRUE,
      colClasses = c("character", rep("numeric", 8))
    )

    quartiles <- summary(lds_seg\$ldscore_SNP)

    lb1 <- which(lds_seg\$ldscore_SNP <= quartiles[2])
    lb2 <- which(lds_seg\$ldscore_SNP > quartiles[2] & lds_seg\$ldscore_SNP <= quartiles[3])
    lb3 <- which(lds_seg\$ldscore_SNP > quartiles[3] & lds_seg\$ldscore_SNP <= quartiles[5])
    lb4 <- which(lds_seg\$ldscore_SNP > quartiles[5])

    write.table(lds_seg\$SNP[lb1], paste(out_prefix, "snp_group1.txt", sep = "_"), row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
    write.table(lds_seg\$SNP[lb2], paste(out_prefix, "snp_group2.txt", sep = "_"), row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
    write.table(lds_seg\$SNP[lb3], paste(out_prefix, "snp_group3.txt", sep = "_"), row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
    write.table(lds_seg\$SNP[lb4], paste(out_prefix, "snp_group4.txt", sep = "_"), row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
    EOF
    """

    stub:
    """
    printf "SNP\tA1\tA2\tFreq\tb\tse\tp\tldscore\n" > ${meta.id}_gcta_ld.score.ld
    printf "stub_snp1\n" > ${meta.id}_snp_group1.txt
    printf "stub_snp2\n" > ${meta.id}_snp_group2.txt
    printf "stub_snp3\n" > ${meta.id}_snp_group3.txt
    printf "stub_snp4\n" > ${meta.id}_snp_group4.txt
    """
}
