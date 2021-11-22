#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR          } from '../../../../modules/untar/main.nf'           addParams( options: [:] )
include { METAPHLAN3     } from '../../../../modules/metaphlan3/main.nf'      addParams( options: [ 'args':'--index mpa_v30_CHOCOPhlAn_201901 --add_viruses --bt2_ps very-sensitive-local' ] )
include { HUMANN_HUMANN  } from '../../../../modules/humann/humann/main.nf'   addParams( options: [:] )
include { HUMANN_REGROUP } from '../../../../modules/humann/regroup/main.nf'  addParams( options: [:] )

workflow test_humann_regroup {

    //METAPHLAN3
    input_metaphlan = [ [ id:'test', single_end:true ], // meta map
              [ file('https://github.com/biobakery/humann/raw/master/examples/demo.fastq.gz', checkIfExists: true) ]
            ]

    db    = channel.fromPath('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/metaphlan_database.tar.gz', type: 'dir', checkIfExists: true)

    UNTAR ( db )
    METAPHLAN3 ( input_metaphlan, UNTAR.out.untar )
    
    //HUMANN3
    input = [ [ [ id:'test', single_end:true ], file('https://github.com/biobakery/humann/raw/master/examples/demo.fastq.gz', checkIfExists: true) ], ]
    ch_input = Channel.from(input)
    nucleotide_db = [ file('http://huttenhower.sph.harvard.edu/humann_data/chocophlan/DEMO_chocophlan.v296_201901b.tar.gz', checkIfExists: true) ]
    protein_db = [ file('http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref90_DEMO_diamond_v201901b.tar.gz', checkIfExists: true) ]
    input_humann = ch_input.join(METAPHLAN3.out.profile)

    HUMANN_HUMANN ( input_humann, nucleotide_db, protein_db )

    HUMANN_REGROUP ( HUMANN_HUMANN.out.genefamilies, 'uniref50_rxn' )
}
