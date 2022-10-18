#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ATLASGENEANNOTATIONMANIPULATION_GTF2FEATUREANNOTATION } from '../../../../../modules/nf-core/atlasgeneannotationmanipulation/gtf2featureannotation/main.nf'

// Make a table of gene annotations

workflow test_atlasgeneannotationmanipulation_gtf2featureannotation {

    gtf_file = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    ATLASGENEANNOTATIONMANIPULATION_GTF2FEATUREANNOTATION ( [[ "id":"homo_sapiens"], gtf_file], [[ "id":"homo_sapiens"], []])
}

// Make a table of transcript annotations with synchonised cDNAs

workflow test_atlasgeneannotationmanipulation_gtf2featureannotation_with_fasta {

    gtf_file = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    transcriptome_fasta = file(params.test_data['homo_sapiens']['genome']['transcriptome_fasta'], checkIfExists: true)

    ATLASGENEANNOTATIONMANIPULATION_GTF2FEATUREANNOTATION ( [[ "id":"homo_sapiens"], gtf_file], [[ "id":"homo_sapiens"], transcriptome_fasta] )
}
