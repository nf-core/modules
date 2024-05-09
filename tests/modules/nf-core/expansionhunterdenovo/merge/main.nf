#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EXPANSIONHUNTERDENOVO_MERGE   } from '../../../../../modules/nf-core/expansionhunterdenovo/merge/main.nf'
include { EXPANSIONHUNTERDENOVO_PROFILE } from '../../../../../modules/nf-core/expansionhunterdenovo/profile/main.nf'

workflow test_expansionhunterdenovo_merge {

    input = Channel.of([
        [ id:'test', type:"case" ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ],
    [
        [ id:'test2', type:"control" ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]
    )

    fasta = [[id:'fasta'],file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]
    fasta_fai = [[id:'fasta_fai'],file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]


    EXPANSIONHUNTERDENOVO_PROFILE (
        input,
        fasta,
        fasta_fai
    )

    manifest = EXPANSIONHUNTERDENOVO_PROFILE.out.str_profile
        .collectFile( name: "manifest.tsv", newLine:true ) { meta, profile ->
                                                ["manifest.tsv", "${meta.id}\t${meta.type}\t${profile}"]
                                            }
        .map({ manifest -> [ [id:"merge"], manifest ]})

    EXPANSIONHUNTERDENOVO_MERGE (
        manifest,
        fasta,
        fasta_fai
    )
}
