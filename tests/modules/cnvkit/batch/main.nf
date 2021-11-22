#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVKIT_BATCH as CNVKIT_HYBRID }      from '../../../../modules/cnvkit/batch/main.nf' addParams( options: [ 'args': '--output-reference reference.cnn' ] )
include { CNVKIT_BATCH as CNVKIT_WGS }         from '../../../../modules/cnvkit/batch/main.nf' addParams( options: [ 'args': '--output-reference reference.cnn --method wgs' ] )
include { CNVKIT_BATCH as CNVKIT_TUMORONLY }   from '../../../../modules/cnvkit/batch/main.nf' addParams( options: [ 'args': '--method wgs' ] )


workflow test_cnvkit_hybrid {
    tumor     = file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    normal    = file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true)

    input     = [ [ id:'test' ], // meta map
                tumor,
                normal
              ]
    fasta     = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    targets   = file(params.test_data['sarscov2']['genome']['baits_bed'], checkIfExists: true)

    CNVKIT_HYBRID ( input, fasta, targets, [] )
}

workflow test_cnvkit_wgs {
    tumor     = file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true) 
    normal    = file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)

    input     = [ [ id:'test'], // meta map
                tumor,
                normal
                ]
    fasta     = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    CNVKIT_WGS ( input, fasta, [], [] )
}


workflow test_cnvkit_cram {
    tumor     = file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)
    normal    = file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)

    input     = [ [ id:'test'], // meta map
                tumor,
                normal
                ]
    fasta     = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    CNVKIT_WGS ( input, fasta, [], [] )
}



workflow test_cnvkit_tumoronly {
    tumor     = file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)

    input     = [ [ id:'test'], // meta map
                tumor,
                [ ]
            ]
    fasta     = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    reference = file(params.test_data['generic']['cnn']['reference'], checkIfExists: true)

    CNVKIT_TUMORONLY ( input, [], [], reference )
}
