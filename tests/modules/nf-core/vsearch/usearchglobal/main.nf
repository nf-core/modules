#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VSEARCH_USEARCHGLOBAL } from '../../../../../modules/nf-core/vsearch/usearchglobal/main.nf'

workflow test_vsearch_usearchglobal {
    
    query = file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true)
    db = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    idcutoff = 0.985
    outoption = "xcfert"  // Nonsense text to check default case.
    columns = "" 
    VSEARCH_USEARCHGLOBAL ( [[id:'test'], query], db, idcutoff, outoption, columns )
}

workflow test_vsearch_usearchglobal_userout {
    
    query = file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true)
    db = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    idcutoff = 0.985
    outoption = "userout"
    columns = "query+target+id" 
    VSEARCH_USEARCHGLOBAL ( [[id:'test'], query], db, idcutoff, outoption, columns )
}
