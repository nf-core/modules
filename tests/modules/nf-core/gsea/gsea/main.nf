#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GSEA_GSEA } from '../../../../../modules/nf-core/gsea/gsea/main.nf'

workflow test_gsea_gsea {
    
    gct_file = file(params.test_data['homo_sapiens']['gene_set_analysis']['gct'], checkIfExists: true)
    cls_file = file(params.test_data['homo_sapiens']['gene_set_analysis']['cls'], checkIfExists: true)
    gmx_file = file(params.test_data['homo_sapiens']['gene_set_analysis']['gmx'], checkIfExists: true)

    input = [ ['id': 'study' ], gct_file, cls_file, gmx_file ]
    contrast = [ 'WT', 'MUT' ] 

    GSEA_GSEA (
        input,
        contrast,
        []
    )
}

workflow test_gsea_gsea_nosets {
    
    gct_file = file(params.test_data['homo_sapiens']['gene_set_analysis']['gct'], checkIfExists: true)
    cls_file = file(params.test_data['homo_sapiens']['gene_set_analysis']['cls'], checkIfExists: true)
    gmx_file = file(params.test_data['homo_sapiens']['gene_set_analysis']['gmx'], checkIfExists: true)

    input = [ ['id': 'study' ], gct_file, cls_file, gmx_file ]
    contrast = [ 'WT', 'MUT' ] 

    GSEA_GSEA (
        input,
        contrast,
        []
    )
}

workflow test_gsea_gsea_zip {
    
    gct_file = file(params.test_data['homo_sapiens']['gene_set_analysis']['gct'], checkIfExists: true)
    cls_file = file(params.test_data['homo_sapiens']['gene_set_analysis']['cls'], checkIfExists: true)
    gmx_file = file(params.test_data['homo_sapiens']['gene_set_analysis']['gmx'], checkIfExists: true)

    input = [ ['id': 'study' ], gct_file, cls_file, gmx_file ]
    contrast = [ 'WT', 'MUT' ] 

    GSEA_GSEA (
        input,
        contrast,
        []
    )
}
