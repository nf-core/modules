include { DEEPVARIANT_RUNDEEPVARIANT      } from '../../../../modules/nf-core/deepvariant/rundeepvariant/main'
include { DEEPVARIANT                     } from '../main'

workflow DEEPVARIANT_WORKFLOW_AND_PROCESS_EQUALITY_TESTER {
    take:
    ch_input   // channel: [ val(meta), path(input), path(index), path(intervals)]
    ch_fasta   // channel: [ val(meta2), path(fasta) ]
    ch_fai     // channel: [ val(meta3), path(fail) ]
    ch_gzi     // channel: [ val(meta4), path(gzi) ]
    ch_par_bed // channel: [ val(meta5), path(par_bed) ]

    main:

    DEEPVARIANT(ch_input, ch_fasta, ch_fai, ch_gzi, ch_par_bed)
    DEEPVARIANT_RUNDEEPVARIANT(ch_input, ch_fasta, ch_fai, ch_gzi, ch_par_bed)

    emit:
    wf_vcf = DEEPVARIANT.out.vcf
    pc_vcf = DEEPVARIANT_RUNDEEPVARIANT.out.vcf
    wf_gvcf = DEEPVARIANT.out.gvcf
    pc_gvcf = DEEPVARIANT_RUNDEEPVARIANT.out.gvcf
}
