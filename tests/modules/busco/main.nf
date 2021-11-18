#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCO } from '../../../modules/busco/main.nf' addParams( options: [args: '--mode genome --auto-lineage'] )

workflow test_busco {
    
   compressed_genome_file = file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)
   input = [ [ id:'test' ], compressed_genome_file ]
   lineage_dataset = file(params.test_data['homo_sapiens']['genome']['chr22_odb10.tar.gz'], checkIfExists: true)

   BUSCO ( input,
            [] )
}
