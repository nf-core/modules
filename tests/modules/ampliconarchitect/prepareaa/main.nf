nextflow.enable.dsl = 2

include { AMPLICONARCHITECT_PREPAREAA } from '../../../../modules/ampliconarchitect/prepareaa/main.nf' addParams( options: [:] )

workflow test_ampliconarchitect_prepareaa {
    input = [ [ id:'test'], //meta map
              file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)
            ]

    BEDTOOLS_BAMTOBED ( input )
}
