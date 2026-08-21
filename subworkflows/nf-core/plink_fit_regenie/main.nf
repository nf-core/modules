//
// Fit REGENIE Step 1 whole-genome regression models from PLINK-format genotype files
//

include { REGENIE_STEP1   } from '../../../modules/nf-core/regenie/step1/main'
include { REGENIE_SPLITL0 } from '../../../modules/nf-core/regenie/splitl0/main'
include { REGENIE_RUNL0   } from '../../../modules/nf-core/regenie/runl0/main'
include { REGENIE_RUNL1   } from '../../../modules/nf-core/regenie/runl1/main'

workflow PLINK_FIT_REGENIE {
    take:
    ch_genotypes // channel: [ val(meta), path(genotype_file), path(variant_file), path(sample_file) ]
    ch_pheno // channel: [ val(meta), path(pheno) ]
    ch_covar // channel: [ val(meta), path(covar) ], use [] when absent
    val_bsize // value: Step 1 block size, use [] to keep module default
    val_step1_mode // value: 'standard' or 'chunked'
    val_n_l0_jobs // value: number of L0 jobs for chunked mode, use [] for standard mode

    main:
    ch_versions = channel.empty()

    if (!['standard', 'chunked'].contains(val_step1_mode)) {
        error("PLINK_FIT_REGENIE: val_step1_mode must be 'standard' or 'chunked', got '${val_step1_mode}'")
    }

    if (val_step1_mode == 'standard') {
        if (!is_empty_value(val_n_l0_jobs)) {
            error("PLINK_FIT_REGENIE: val_n_l0_jobs must be empty when val_step1_mode is 'standard'")
        }

        REGENIE_STEP1(
            ch_genotypes,
            ch_pheno,
            ch_covar,
            val_bsize,
        )

        ch_predictions = REGENIE_STEP1.out.predictions
        ch_loco = REGENIE_STEP1.out.loco
        ch_log = REGENIE_STEP1.out.log
    }
    else {
        n_l0_jobs = require_positive_integer(
            val_n_l0_jobs,
            "PLINK_FIT_REGENIE: val_n_l0_jobs must be a positive integer when val_step1_mode is 'chunked'",
        )

        REGENIE_SPLITL0(
            ch_genotypes,
            ch_pheno,
            ch_covar,
            val_bsize,
            n_l0_jobs,
        )

        ch_l0_jobs = REGENIE_SPLITL0.out.master
            .join(REGENIE_SPLITL0.out.snplists, by: [0])
            .flatMap { meta, master, snplists ->
                (1..n_l0_jobs).collect { job_number ->
                    [meta, master, snplists, job_number]
                }
            }

        ch_runl0_inputs = ch_genotypes
            .combine(ch_l0_jobs, by: [0])
            .combine(ch_pheno, by: [0])
            .combine(ch_covar, by: [0])

        REGENIE_RUNL0(
            ch_runl0_inputs.map { values ->
                [values[0], values[1], values[2], values[3]]
            },
            ch_runl0_inputs.map { values ->
                [values[0], values[4], values[5], values[6]]
            },
            ch_runl0_inputs.map { values ->
                [values[0], values[7]]
            },
            ch_runl0_inputs.map { values ->
                [values[0], values[8]]
            },
            val_bsize,
        )

        ch_l0_predictions = REGENIE_RUNL0.out.l0_predictions
            .groupTuple(by: [0])
            .map { meta, l0_predictions ->
                [meta, l0_predictions.flatten()]
            }

        ch_runl1_split_inputs = REGENIE_SPLITL0.out.master
            .join(REGENIE_SPLITL0.out.snplists, by: [0])
            .join(ch_l0_predictions, by: [0])

        ch_runl1_inputs = ch_genotypes
            .combine(ch_runl1_split_inputs, by: [0])
            .combine(ch_pheno, by: [0])
            .combine(ch_covar, by: [0])

        REGENIE_RUNL1(
            ch_runl1_inputs.map { values ->
                [values[0], values[1], values[2], values[3]]
            },
            ch_runl1_inputs.map { values ->
                [values[0], values[4], values[5], values[6]]
            },
            ch_runl1_inputs.map { values ->
                [values[0], values[7]]
            },
            ch_runl1_inputs.map { values ->
                [values[0], values[8]]
            },
            val_bsize,
        )

        ch_predictions = REGENIE_RUNL1.out.predictions
        ch_loco = REGENIE_RUNL1.out.loco
        ch_log = REGENIE_SPLITL0.out.log.mix(REGENIE_RUNL0.out.log, REGENIE_RUNL1.out.log)
    }

    emit:
    predictions = ch_predictions
    loco        = ch_loco
    logs        = ch_log
    versions    = ch_versions
}

def is_empty_value(value) {
    value == null || value == [] || value == ''
}

def require_positive_integer(value, message) {
    try {
        def integer_value = value as Integer
        if (integer_value < 1) {
            error(message)
        }
        return integer_value
    }
    catch (_ignored) {
        error(message)
    }
}
