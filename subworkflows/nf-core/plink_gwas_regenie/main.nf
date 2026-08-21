//
// Run REGENIE Step 1 fitting and Step 2 association testing from PLINK-format genotype files
//

include { PLINK_FIT_REGENIE } from '../plink_fit_regenie/main'
include { REGENIE_STEP2     } from '../../../modules/nf-core/regenie/step2/main'

workflow PLINK_GWAS_REGENIE {
    take:
    ch_step1_genotypes // channel: [ val(meta), path(genotype_file), path(variant_file), path(sample_file) ]
    ch_step2_genotypes // channel: [ val(meta), path(genotype_file), path(variant_file), path(sample_file) ]
    ch_pheno // channel: [ val(meta), path(pheno) ]
    ch_covar // channel: [ val(meta), path(covar) ], use [] when absent
    val_step1_bsize // value: Step 1 block size, use [] to keep module default
    val_step2_bsize // value: Step 2 block size
    val_step1_mode // value: 'standard' or 'chunked'
    val_n_l0_jobs // value: number of L0 jobs for chunked mode, use [] for standard mode

    main:
    ch_versions = channel.empty()

    PLINK_FIT_REGENIE(
        ch_step1_genotypes,
        ch_pheno,
        ch_covar,
        val_step1_bsize,
        val_step1_mode,
        val_n_l0_jobs,
    )

    ch_step1_predictions = PLINK_FIT_REGENIE.out.predictions.map { meta, predictions ->
        [meta_id(meta, 'ch_step1_genotypes'), meta, predictions]
    }

    ch_step1_loco = PLINK_FIT_REGENIE.out.loco.map { meta, loco ->
        [meta_id(meta, 'ch_step1_genotypes'), meta, loco]
    }

    ch_step1_pred_loco = ch_step1_predictions
        .join(ch_step1_loco, by: [0])
        .ifEmpty {
            error("PLINK_GWAS_REGENIE: no matching Step 1 predictions and LOCO files were available for REGENIE_STEP2")
        }
        .map { values ->
            [values[0], values[1], values[2], values[4]]
        }

    ch_keyed_step2_genotypes = ch_step2_genotypes.map { meta, genotype_file, variant_file, sample_file ->
        [meta_id(meta, 'ch_step2_genotypes'), meta, genotype_file, variant_file, sample_file]
    }

    ch_keyed_pheno = ch_pheno.map { meta, pheno ->
        [meta_id(meta, 'ch_pheno'), meta, pheno]
    }

    ch_keyed_covar = ch_covar.map { meta, covar ->
        [meta_id(meta, 'ch_covar'), meta, covar]
    }

    ch_step2_inputs = ch_keyed_step2_genotypes
        .combine(ch_step1_pred_loco, by: [0])
        .ifEmpty {
            error("PLINK_GWAS_REGENIE: no Step 2 genotype records matched Step 1 predictions by meta.id")
        }
        .combine(ch_keyed_pheno, by: [0])
        .ifEmpty {
            error("PLINK_GWAS_REGENIE: no Step 2 records matched phenotype records by meta.id")
        }
        .combine(ch_keyed_covar, by: [0])
        .ifEmpty {
            error("PLINK_GWAS_REGENIE: no Step 2 records matched covariate records by meta.id")
        }

    REGENIE_STEP2(
        ch_step2_inputs.map { values ->
            [values[1], values[2], values[3], values[4]]
        },
        ch_step2_inputs.map { values ->
            [values[5], values[6], values[7]]
        },
        ch_step2_inputs.map { values ->
            [values[8], values[9]]
        },
        ch_step2_inputs.map { values ->
            [values[10], values[11]]
        },
        val_step2_bsize,
    )

    emit:
    results     = REGENIE_STEP2.out.results
    logs        = PLINK_FIT_REGENIE.out.logs.mix(REGENIE_STEP2.out.log)
    predictions = PLINK_FIT_REGENIE.out.predictions
    loco        = PLINK_FIT_REGENIE.out.loco
    versions    = ch_versions
}

def meta_id(meta, channel_name) {
    if (!meta?.id) {
        error("PLINK_GWAS_REGENIE: ${channel_name} records must include meta.id")
    }
    meta.id
}
