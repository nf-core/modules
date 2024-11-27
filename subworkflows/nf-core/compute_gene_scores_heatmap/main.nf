//
// Compute gene scores and plot heatmap
//

//
// MODULES
//
include { COMPUTE_GENE_SCORES     } from '../../../modules/nf-core/compute_gene_scores'
include { CREATE_GENE_HEATMAP     } from '../../../modules/nf-core/create_gene_heatmap'

workflow COMPUTE_GENE_SCORES_HEATMAP {
    take:
    normalized_counts          // channel: [ meta, normalized_counts.tsv ]
    annotated_endo_data        // channel: [ meta, annotated_endo_data.tsv ]
    ch_gene_score_yaml         // channel: [ file(gene_score_yaml.yaml) ]
    ch_heatmap_genes_to_filter // channel: [ file(heatmap_genes_to_filter.yaml) ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Compute gene scores for supplied YAML gene score file
    //
    COMPUTE_GENE_SCORES(
        normalized_counts,
        ch_gene_score_yaml
    )
    ch_versions      = ch_versions.mix(COMPUTE_GENE_SCORES.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(COMPUTE_GENE_SCORES.out.scores_for_mqc.map{ meta, file-> file }.collect())

    //
    // MODULE: Compute gene-count heatmap for MultiQC report based on annotated (ENDO) counts
    //
    if(!params.skip_heatmap){
        ch_create_gene_heatmap_input = annotated_endo_data.join(normalized_counts)
        CREATE_GENE_HEATMAP (
            ch_create_gene_heatmap_input,
            ch_heatmap_genes_to_filter.toList()
        )
        ch_versions       = ch_versions.mix(CREATE_GENE_HEATMAP.out.versions)
        ch_multiqc_files  = ch_multiqc_files.mix(CREATE_GENE_HEATMAP.out.gene_heatmap.map{ meta, file-> file }.collect())
    }

    emit:
    versions                = ch_versions      // channel: [ versions.yml ]
    multiqc_files           = ch_multiqc_files // channel: [*mqc.txt, *gene_heatmap_mqc.png ]
}
