//
// Infer protein activity from network interactions of a set of
//     gene regulators from gene expression data
//

include { ARACNE3 } from '../../../modules/nf-core/aracne3/main'
include { VIPER } from '../../../modules/nf-core/viper/main'

workflow EXPRESSION_PROTEINACTIVITY {

    take:
    ch_regulators        // channel: /path/to/regulators_genes.txt
    ch_expression_matrix // channel: /path/to/expression_quant.sf
    ch_network           // channel: /path/to/regulon_network.txt
    ch_sample_expression // channel: [ val(meta), [ reads ] ]

    main:

    ch_versions = Channel.empty()

    // Run ARACNE
    if (!ch_network) {
        ch_network = ARACNE3 ( ch_expression_matrix, ch_regulators ).consensus_network
        ch_versions = ch_versions.mix(ARACNE3.out.versions)
    }

    // Run VIPER
    VIPER ( ch_network, ch_sample_expression )
    ch_versions = ch_versions.mix(VIPER.out.versions)

    emit:
    aracne_network = ARACNE3.out.consensus_network      // channel: [ tsv ]
    aracne_subnets = ARACNE3.out.subnets                // channel: [ dir ]
    viper_results  = VIPER.out.viper_results            // channel: [ tsv ]
    versions       = ch_versions                        // channel: [ versions.yml ]
}


