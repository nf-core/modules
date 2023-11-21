#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MUDSKIPPER_BULK  } from '../../../../../modules/nf-core/mudskipper/bulk/main.nf'
include { MUDSKIPPER_INDEX } from '../../../../../modules/nf-core/mudskipper/index/main.nf'
include { SAMTOOLS_SORT    } from '../../../../../modules/nf-core/samtools/sort/main.nf'

workflow test_mudskipper_bulk_gtf {

    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    
    input = [
        [id:'test', single_end:false ], // meta map
        file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test.rna.paired_end.sorted.bam", checkIfExists: true)
    ]
    
    SAMTOOLS_SORT ( input ) 

    MUDSKIPPER_BULK ( SAMTOOLS_SORT.out.bam, [], gtf, false )
}

workflow test_mudskipper_bulk_index {

    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    
    input = [
        [id:'test', single_end:false ], // meta map
        file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test.rna.paired_end.sorted.bam", checkIfExists: true)
    ]

    SAMTOOLS_SORT ( input )

    MUDSKIPPER_INDEX ( gtf )
    MUDSKIPPER_BULK ( SAMTOOLS_SORT.out.bam, MUDSKIPPER_INDEX.out.index, [], false )

}
