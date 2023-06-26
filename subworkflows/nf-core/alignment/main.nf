// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { BWA_MEM      } from '../../../modules/nf-core/bwa/mem/main'
include { BWAMEM2_MEM      } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_INDEX      } from '../../../modules/nf-core/bwamem2/index/main'
include { PICARD_ADDORREPLACEREADGROUPS      } from '../../../modules/nf-core/picard/addorreplacereadgroups/main'

workflow ALIGNMENT {

    take:
    // TODO nf-core: edit input (take) channels
    fastqs // channel: [ val(meta), [ bam ] ]
    reference
    bwa

    main:

    versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow
    if(bwa == 1) {
        BWA_MEM ( fastqs, reference, true ).bam.map {
            meta, bam ->
                new_id = 'aligned_bam'
                [[id: new_id], bam ]
        }.set {aligned_bam}
        versions = versions.mix(BWA_MEM.out.versions.first())
    } 
    if(bwa == 2) {
        // Index Reference fasta
        fasta_file = reference[1].list().find{it=~/.fasta$/}
        fasta_path = [
        [id:'reference_fasta'], file(reference[1].resolve(fasta_file))
        ]
        BWAMEM2_INDEX (fasta_path)
        versions = versions.mix(BWAMEM2_INDEX.out.versions.first())

        // BWA MEM2 
        BWAMEM2_MEM ( fastqs, BWAMEM2_INDEX.out.index, true ).bam.map {
            meta, bam ->
                new_id = 'aligned_bam'
                [[id: new_id], bam ]
        }.set {aligned_bam}
        versions = versions.mix(BWAMEM2_MEM.out.versions.first())
    }


    PICARD_ADDORREPLACEREADGROUPS(aligned_bam).bam.map {
        meta, bam ->
            new_id = 'grouped_aligned_bam'
            [[id: new_id], bam ]
    }.set {grouped_bam}
    versions = versions.mix(PICARD_ADDORREPLACEREADGROUPS.out.versions.first())
    emit:
    // TODO nf-core: edit emitted channels
    bam      = PICARD_ADDORREPLACEREADGROUPS.out.bam           // channel: [ val(meta), [ bam ] ]


    versions = versions                     // channel: [ versions.yml ]
}

