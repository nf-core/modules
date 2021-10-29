#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

test_filters = "--attribute=gene_biotype:protein_coding \
                    --attribute=gene_biotype:lincRNA \
                    --attribute=gene_biotype:antisense \
                    --attribute=gene_biotype:IG_LV_gene \
                    --attribute=gene_biotype:IG_V_gene \
                    --attribute=gene_biotype:IG_V_pseudogene \
                    --attribute=gene_biotype:IG_D_gene \
                    --attribute=gene_biotype:IG_J_gene \
                    --attribute=gene_biotype:IG_J_pseudogene \
                    --attribute=gene_biotype:IG_C_gene \
                    --attribute=gene_biotype:IG_C_pseudogene \
                    --attribute=gene_biotype:TR_V_gene \
                    --attribute=gene_biotype:TR_V_pseudogene \
                    --attribute=gene_biotype:TR_D_gene \
                    --attribute=gene_biotype:TR_J_gene \
                    --attribute=gene_biotype:TR_J_pseudogene \
                    --attribute=gene_biotype:TR_C_gene"

include { CELLRANGER_MKGTF } from '../../../../modules/cellranger/mkgtf/main.nf' addParams( options: [args: test_filters] )

workflow test_cellranger_mkgtf {
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    CELLRANGER_MKGTF ( gtf )
}
