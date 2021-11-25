#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CREATE_UMI_CONSENSUS } from '../../../../subworkflows/nf-core/fgbio_create_umi_consensus/main' addParams( [:]  )
params.read_structure = "+T 12M11S+T"

workflow test_fgbio_create_umi_consensus {
    reads          = [ [ id:'test', single_end:false ], // meta map
                       [ file(params.test_data['homo_sapiens']['illumina']['test_umi_1_fastq_gz'], checkIfExists: true),
                        file(params.test_data['homo_sapiens']['illumina']['test_umi_2_fastq_gz'],  checkIfExists: true) ]
                     ]
    fasta          = file(params.test_data['homo_sapiens']['genome']['genome_fasta'],              checkIfExists: true)

    read_structure =

    CREATE_UMI_CONSENSUS( reads, fasta, params.read_structure)
}
