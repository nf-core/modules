#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR        } from '../../../../../modules/nf-core/untar/main.nf'
include { AFFY_JUSTRMA } from '../../../../../modules/nf-core/affy/justrma/main.nf'
include { LIMMA_DIFFERENTIAL } from '../../../../../modules/nf-core/limma/differential/main.nf'

workflow test_limma_differential {

    // Make a microarray intensities matrix using affy/justrma

    samples = file(params.test_data['homo_sapiens']['genome']['affy_array_samplesheet'], checkIfExists: true)
    cel_archive = file(params.test_data['homo_sapiens']['genome']['affy_array_celfiles_tar'], checkIfExists: true)
    
    def meta = [ id:'test' ]
    ch_samplesheet = Channel.of([ meta, samples ])
    ch_celfiles_archive = Channel.of([ meta, cel_archive ])

    UNTAR ( ch_celfiles_archive )

    ch_input = ch_samplesheet.join(UNTAR.out.untar)

    AFFY_JUSTRMA ( 
        ch_input,
        [[],[]] 
    )

    // Now do the actual differential analysis

    ch_intensities = ch_samplesheet
        .join(AFFY_JUSTRMA.out.expression)

    ch_contrasts = Channel.of(['id': 'diagnosis_normal_uremia', 'variable': 'diagnosis', 'reference': 'normal', 'target': 'uremia'])
        .map{
            tuple(it, it.variable, it.reference, it.target)
        }

    LIMMA_DIFFERENTIAL (
        ch_contrasts,
        ch_intensities
    )
}

workflow test_limma_differential_subset_to_contrast {

    // Make a microarray intensities matrix using affy/justrma

    samples = file(params.test_data['homo_sapiens']['genome']['affy_array_samplesheet'], checkIfExists: true)
    cel_archive = file(params.test_data['homo_sapiens']['genome']['affy_array_celfiles_tar'], checkIfExists: true)
    
    def meta = [ id:'test' ]
    ch_samplesheet = Channel.of([ meta, samples ])
    ch_celfiles_archive = Channel.of([ meta, cel_archive ])

    UNTAR ( ch_celfiles_archive )

    ch_input = ch_samplesheet.join(UNTAR.out.untar)

    AFFY_JUSTRMA ( 
        ch_input,
        [[],[]] 
    )

    // Now do the actual differential analysis

    ch_intensities = ch_samplesheet
        .join(AFFY_JUSTRMA.out.expression)

    ch_contrasts = Channel.of(['id': 'diagnosis_normal_uremia', 'variable': 'diagnosis', 'reference': 'normal', 'target': 'uremia'])
        .map{
            tuple(it, it.variable, it.reference, it.target)
        }

    LIMMA_DIFFERENTIAL (
        ch_contrasts,
        ch_intensities
    )
}


workflow test_limma_differential_exclude_samples {

    // Make a microarray intensities matrix using affy/justrma

    samples = file(params.test_data['homo_sapiens']['genome']['affy_array_samplesheet'], checkIfExists: true)
    cel_archive = file(params.test_data['homo_sapiens']['genome']['affy_array_celfiles_tar'], checkIfExists: true)
    
    def meta = [ id:'test' ]
    ch_samplesheet = Channel.of([ meta, samples ])
    ch_celfiles_archive = Channel.of([ meta, cel_archive ])

    UNTAR ( ch_celfiles_archive )

    ch_input = ch_samplesheet.join(UNTAR.out.untar)

    AFFY_JUSTRMA ( 
        ch_input,
        [[],[]] 
    )

    // Now do the actual differential analysis

    ch_intensities = ch_samplesheet
        .join(AFFY_JUSTRMA.out.expression)

    ch_contrasts = Channel.of(['id': 'diagnosis_normal_uremia', 'variable': 'diagnosis', 'reference': 'normal', 'target': 'uremia'])
        .map{
            tuple(it, it.variable, it.reference, it.target)
        }

    LIMMA_DIFFERENTIAL (
        ch_contrasts,
        ch_intensities
    )
}
