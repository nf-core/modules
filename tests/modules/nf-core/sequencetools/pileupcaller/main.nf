#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQUENCETOOLS_PILEUPCALLER as SEQUENCETOOLS_PILEUPCALLER_FREQSUM } from '../../../../../modules/nf-core/sequencetools/pileupcaller/main.nf'
include { SEQUENCETOOLS_PILEUPCALLER as SEQUENCETOOLS_PILEUPCALLER_EIG } from '../../../../../modules/nf-core/sequencetools/pileupcaller/main.nf'
include { SEQUENCETOOLS_PILEUPCALLER as SEQUENCETOOLS_PILEUPCALLER_PLINK } from '../../../../../modules/nf-core/sequencetools/pileupcaller/main.nf'

workflow test_sequencetools_pileupcaller {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_mpileup'], checkIfExists: true)
            ]
    snpfile = file(params.test_data['homo_sapiens']['genome']['genome_21_eigenstrat_snp'], checkIfExists: true)

    SEQUENCETOOLS_PILEUPCALLER_FREQSUM ( input, snpfile, [] )
    SEQUENCETOOLS_PILEUPCALLER_EIG ( input, snpfile, [] )
    SEQUENCETOOLS_PILEUPCALLER_PLINK ( input, snpfile, [] )
}
