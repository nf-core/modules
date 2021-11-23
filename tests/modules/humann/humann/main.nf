#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR          } from '../../../../modules/untar/main.nf'           addParams( options: [:] )
include { SAMTOOLS_VIEW  } from '../../../../modules/samtools/view/main.nf'   addParams( options: ['suffix': '.sam'] )
include { METAPHLAN3     } from '../../../../modules/metaphlan3/main.nf'      addParams( options: [ 'args':'--index mpa_v30_CHOCOPhlAn_201901 --bt2_ps very-sensitive-local' ] )
include { HUMANN_HUMANN  } from '../../../../modules/humann/humann/main.nf'   addParams( options: [:] )

process RENAME {
    input:
    path file_in
    val  file_out

    output:
    path "${file_out}*" , emit: file_out

    script:
    """
    mv $file_in $file_out
    """
}

workflow test_humann_humann {
    //Rename input file because only 'demo.fastq.gz' is accepted with DEMO (i.e. small) databases
    RENAME ( file(params.test_data['bacteroides_fragilis']['illumina']['test1_1_fastq_gz'], checkIfExists: true), 'demo.fastq.gz' )
    
    //Create meta map
    input = [ [ id:'test', single_end:true ], ]
    ch_input = Channel.from(input)

    //METAPHLAN3
    input_metaphlan = ch_input.combine(RENAME.out.file_out)
    db    = channel.fromPath('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/metaphlan_database.tar.gz', type: 'dir', checkIfExists: true)

    UNTAR ( db )
    METAPHLAN3 ( input_metaphlan, UNTAR.out.untar )
    
    //HUMANN3
    input_humann = ch_input.combine(RENAME.out.file_out)
    input_humann = input_humann.join(METAPHLAN3.out.profile)
    nucleotide_db = [ file('http://huttenhower.sph.harvard.edu/humann_data/chocophlan/DEMO_chocophlan.v296_201901b.tar.gz', checkIfExists: true) ]
    protein_db = [ file('http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref90_DEMO_diamond_v201901b.tar.gz', checkIfExists: true) ]

    HUMANN_HUMANN ( input_humann, nucleotide_db, protein_db )
}
