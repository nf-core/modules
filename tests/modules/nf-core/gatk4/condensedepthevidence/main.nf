#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_CONDENSEDEPTHEVIDENCE } from '../../../../../modules/nf-core/gatk4/condensedepthevidence/main.nf'

workflow test_gatk4_condensdepthevidence {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/condensedepthevidence/testN.rd.txt.gz", checkIfExists: true),
        file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/condensedepthevidence/testN.rd.txt.gz.tbi", checkIfExists: true)
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_CONDENSEDEPTHEVIDENCE (
        input,
        fasta,
        fasta_fai,
        dict
    )
}
