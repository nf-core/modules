/*
 * Runs of homozygosity (ROH) analysis, same steps as the Goat ROH pipeline
 * (main.py of the original Goat ROH analysis):
 *   1. PLINK --homozyg                       -> .hom, .hom.indiv, .hom.summary, .log
 *   2. F_ROH per FID                         -> F_ROH_estimates.csv
 *   3. IID -> breed metadata from .fam       -> SheepAnalysis_ROH_metadata.csv
 *   4. ROH counts per length bin             -> ROH_bins_per_individual.csv, ROH_bins_per_breed.csv
 *   5. ROH total length per length bin       -> ROH_total_length_per_individual.csv, ROH_total_length_per_breed.csv
 *   6. Count-vs-length and composition plots -> F_ROH_per_individual.csv, plots
 *
 * Defaults are PLINK's own --homozyg defaults, as in the original analysis. They
 * are not tuned for a 60k goat array; pass explicit --homozyg-* values in ext.args
 * to change them, and they will be recorded in the log.
 *
 * Genome length: both steps default to 2,500 Mb (2,500,000 kb), matching the recorded
 * Goat results. Change it per run with ext.genome_kb (step 2) and ext.genome_mb (step 6),
 * e.g. in a config: process { withName: 'PLINK_HOMOZYG' { ext.genome_mb = 2452 } }
 * Both are kept as in the original analysis; override with ext.genome_kb / ext.genome_mb.
 *
 * Note: F_ROH from step 6 sums only the five length bins (>= 1 Mb). ROH shorter
 * than 1 Mb are not included, as in the original analysis.
 *
 * Container: this process needs PLINK and Python together. The PLINK-only image on the
 * line below does not contain Python, so for the Goat validation the process was run
 * with the combined image built from environment.yml:
 *   validation/container/roh_homozyg.sif (set in validation/roh_full/nextflow.config).
 * Publish that combined image (mulled container) before upstream submission and
 * update the container line.
 *
 * Changes from the original scripts (logic otherwise unchanged):
 *   - 6_visualize_roh.py: a bin with no ROH in the dataset is added as zero, so
 *     datasets without, for example, 1-2 Mb ROH no longer fail at step 6.
 *   - 3_create_metadata.py: output name is an argument (default unchanged).
 *   - All scripts: delim_whitespace=True replaced by sep=r"\\s+" (same behaviour; pandas 3 removed the old argument).
 *
 * Validation (Goat SaudiGoatArdi, repositioned input, 144 animals, 59,727 SNPs):
 *   - .hom, .hom.indiv, .hom.summary identical to the recorded run (188 ROH).
 *   - F_ROH, ROH bin counts and breed tables identical.
 *   - Step 6 default changed from the script's 2,452 Mb to 2,500 Mb, so F_ROH matches
 *     the recorded 6.F_ROH_per_individual.csv from the original analysis.
 *     that file.
 *   - ROH_total_length_per_individual.csv: one row of the recorded file disagrees
 *     with the recorded bin counts of the same run; this module agrees with the counts.
 */
process PLINK_HOMOZYG {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Combined PLINK 1.90b6.21 + Python image (built with mulled-build from this module's environment.yml).
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/othmanaljurayyad/plink-homozyg:1.90b6.21' :
        'ghcr.io/othmanaljurayyad/plink-homozyg:1.90b6.21' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("*.hom")                       , emit: hom
    tuple val(meta), path("*.hom.indiv")                 , emit: hom_indiv
    tuple val(meta), path("*.hom.summary")               , emit: hom_summary
    tuple val(meta), path("*.log")                       , emit: log
    tuple val(meta), path("F_ROH_estimates.csv")         , emit: froh_breed
    tuple val(meta), path("F_ROH_per_individual.csv")    , emit: froh_individual
    tuple val(meta), path("*_metadata.csv")              , emit: metadata
    tuple val(meta), path("ROH_bins_per_individual.csv") , emit: bins_individual
    tuple val(meta), path("ROH_bins_per_breed.csv")      , emit: bins_breed
    tuple val(meta), path("ROH_total_length_per_individual.csv"), emit: length_individual
    tuple val(meta), path("ROH_total_length_per_breed.csv"),      emit: length_breed
    tuple val(meta), path("*.png")                       , emit: plots
    tuple val("${task.process}"), val('plink'), eval("plink --version 2>&1 | sed 's/^PLINK v//;s/ .*//'"), emit: versions_plink, topic: versions
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/^Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def genome_kb   = task.ext.genome_kb ?: 2500000
    def genome_mb   = task.ext.genome_mb ?: 2500
    """
    plink \\
        --bfile ${bed.baseName} \\
        --homozyg \\
        $args \\
        --threads $task.cpus \\
        --out ${prefix}_ROH

    python3 ${moduleDir}/bin/2_calculate_froh.py  --roh ${prefix}_ROH.hom --genome_kb ${genome_kb} --out .
    python3 ${moduleDir}/bin/3_create_metadata.py --fam ${fam} --out . --name ${prefix}_ROH_metadata.csv
    python3 ${moduleDir}/bin/4_bin_roh_counts.py  --roh ${prefix}_ROH.hom --meta ${prefix}_ROH_metadata.csv --out .
    python3 ${moduleDir}/bin/5_bin_roh_lengths.py --roh ${prefix}_ROH.hom --meta ${prefix}_ROH_metadata.csv --out .
    python3 ${moduleDir}/bin/6_visualize_roh.py   --roh_counts ROH_bins_per_individual.csv --roh_lengths ROH_total_length_per_individual.csv --out . --genome_mb ${genome_mb}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_ROH.hom ${prefix}_ROH.hom.indiv ${prefix}_ROH.hom.summary ${prefix}_ROH.log
    touch F_ROH_estimates.csv F_ROH_per_individual.csv ${prefix}_ROH_metadata.csv
    touch ROH_bins_per_individual.csv ROH_bins_per_breed.csv
    touch ROH_total_length_per_individual.csv ROH_total_length_per_breed.csv
    touch plot.png
    """
}
