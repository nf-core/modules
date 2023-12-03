#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {VGAN_HAPLOCART} from '../../../../../modules/nf-core/vgan/haplocart/main.nf'

process STUB_HCFILES {

    output:
        path "hcfiles", emit: hcfiles

    script:
        """
            mkdir hcfiles
            touch hcfiles/graph.gg
            touch hcfiles/graph.og
            touch hcfiles/graph.xg
            touch hcfiles/graph.giraffe.gbz
            touch hcfiles/graph.dist
            touch hcfiles/graph_paths
            touch hcfiles/path_supports
            touch hcfiles/graph.gbwt
            touch hcfiles/children.txt
            touch hcfiles/parents.txt
            touch hcfiles/parsed_pangenome_mapping
            touch hcfiles/mappability.tsv
            touch hcfiles/k17_w18.min
            touch hcfiles/k31_w11.min
    """
}


workflow test_vgan_haplocart_single_end {

    //ch_vg_db = STUB_HCFILES().hcfiles
    ch_vg_db = file("/Users/andrades/vgan/hcfiles")

    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['rCRS_reads'], checkIfExists: true) ]
            ]

    VGAN_HAPLOCART(input, ch_vg_db, [])
}

workflow test_vgan_haplocart_paired_end_separate {

    //ch_vg_db = STUB_HCFILES().hcfiles
    ch_vg_db = file("/Users/andrades/vgan/hcfiles")


    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['rCRS_reads'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['rCRS_reads'], checkIfExists: true) ]
            ]

    VGAN_HAPLOCART(input, ch_vg_db, [])
}

workflow test_vgan_haplocart_paired_end_interleaved {

    //ch_vg_db = STUB_HCFILES().hcfiles
    ch_vg_db = file("/Users/andrades/vgan/hcfiles")

    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['rCRS_reads'], checkIfExists: true) ]
            ]

    interleaved =  true

    VGAN_HAPLOCART(input, STUB_HCFILES.out.hcfiles, interleaved)

}
