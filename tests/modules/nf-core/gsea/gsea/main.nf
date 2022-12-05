#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GSEA_GSEA } from '../../../../../modules/nf-core/gsea/gsea/main.nf'

workflow test_gsea_gsea {
    
    gct_file = file(params.test_data['homo_sapiens']['gene_set_analysis']['gct'], checkIfExists: true)
    cls_file = file(params.test_data['homo_sapiens']['gene_set_analysis']['cls'], checkIfExists: true)
    gmx_file = file(params.test_data['homo_sapiens']['gene_set_analysis']['gmx'], checkIfExists: true)

    input = [ ['id': 'study', 'variable': 'treatment', 'reference': 'WT', 'target': 'MUT', 'gene_set_db': 'test_gene_set'], gct_file, cls_file, gmx_file ]

    GSEA_GSEA (
        input 
    )
}

workflow test_gsea_gsea_nosets {
    
    gct_file = file(params.test_data['homo_sapiens']['gene_set_analysis']['gct'], checkIfExists: true)
    cls_file = file(params.test_data['homo_sapiens']['gene_set_analysis']['cls'], checkIfExists: true)
    gmx_file = file(params.test_data['homo_sapiens']['gene_set_analysis']['gmx'], checkIfExists: true)

    input = [ ['id': 'study', 'variable': 'treatment', 'reference': 'WT', 'target': 'MUT', 'gene_set_db': 'test_gene_set'], gct_file, cls_file, gmx_file ]

    GSEA_GSEA (
        input 
    )
}

workflow test_gsea_gsea_zip {
    
    gct_file = file(params.test_data['homo_sapiens']['gene_set_analysis']['gct'], checkIfExists: true)
    cls_file = file(params.test_data['homo_sapiens']['gene_set_analysis']['cls'], checkIfExists: true)
    gmx_file = file(params.test_data['homo_sapiens']['gene_set_analysis']['gmx'], checkIfExists: true)

    input = [ ['id': 'study', 'variable': 'treatment', 'reference': 'WT', 'target': 'MUT', 'gene_set_db': 'test_gene_set'], gct_file, cls_file, gmx_file ]

    GSEA_GSEA (
        input 
    )
}
