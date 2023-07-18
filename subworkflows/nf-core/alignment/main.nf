// import
// bwa2 for extra alignment option
include { BWA_MEM      } from '../../../modules/nf-core/bwa/mem/main'
include { BWA_INDEX      } from '../../../modules/nf-core/bwa/index/main'
include { BWAMEM2_MEM      } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_INDEX      } from '../../../modules/nf-core/bwamem2/index/main'
include { PICARD_ADDORREPLACEREADGROUPS      } from '../../../modules/nf-core/picard/addorreplacereadgroups/main'

workflow ALIGNMENT {

    take:
    fastqs // channel: [ val(meta), [ bam ] ]
    reference
    bwa

    main:

    versions = Channel.empty()

    // switch statement to determine which bwa to use, this is a passed parameter
    switch(bwa){
        case 1:
            // INDEX
            BWA_INDEX ( reference )
            versions = versions.mix(BWA_INDEX.out.versions.first())
            // MEM
            BWA_MEM ( fastqs, BWA_INDEX.out.index, true ).bam.map {
                meta, bam ->
                    new_id = 'aligned_bam'
                    [[id: new_id], bam ]
            }.set {aligned_bam}
            versions = versions.mix(BWA_MEM.out.versions.first())
            break
        case 2:
            // INDEX
            BWAMEM2_INDEX (reference)
            versions = versions.mix(BWAMEM2_INDEX.out.versions.first())
            // BWA MEM2
            BWAMEM2_MEM ( fastqs, BWAMEM2_INDEX.out.index, true ).bam.map {
                meta, bam ->
                    new_id = 'aligned_bam'
                    [[id: new_id], bam ]
            }.set {aligned_bam}
            versions = versions.mix(BWAMEM2_MEM.out.versions.first())
            break
        default:
            throw new Exception("The argument bwa must be either 1 or 2, not ${bwa}.")
    }


    PICARD_ADDORREPLACEREADGROUPS(aligned_bam).bam.map {
        meta, bam ->
            new_id = 'grouped_aligned_bam'
            [[id: new_id], bam ]
    }.set {grouped_bam}
    versions = versions.mix(PICARD_ADDORREPLACEREADGROUPS.out.versions.first())

    // final output
    emit:

    bam      = PICARD_ADDORREPLACEREADGROUPS.out.bam           // channel: [ val(meta), [ bam ] ]


    versions = versions                     // channel: [ versions.yml ]
}

