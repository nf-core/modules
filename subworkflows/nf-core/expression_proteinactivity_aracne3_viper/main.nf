//
// Infer protein activity from network interactions of a set of
//     gene regulators from gene expression data
//

include { ARACNE3 } from '../../../modules/nf-core/aracne3/main'
include { VIPER } from '../../../modules/nf-core/viper/main'

workflow EXPRESSION_PROTEINACTIVITY {

    take:
    ch_expression_matrix // channel: [ val(meta), /path/to/expression.tsv ]
    ch_regulators        // channel: /path/to/regulators.txt
    ch_network           // channel: [ val(meta), /path/to/network.tsv ]
    ch_expression_sample // channel: [ val(meta), /path/to/expression.tsv ]

    main:

    ch_versions = Channel.empty()

    // Run ARACNE
    if (ch_network.empty) {
        ch_network = ARACNE3 ( ch_expression_matrix, ch_regulators ).consensus_network
        ch_versions = ch_versions.mix(ARACNE3.out.versions)
    }

    // Run VIPER
    if (ch_expression_sample.empty) {
        ch_expression_sample = ch_expression_matrix
    }
    VIPER ( ch_network, ch_expression_sample )
    ch_versions = ch_versions.mix(VIPER.out.versions)

    emit:
    aracne_network = ARACNE3.out.consensus_network      // channel: [ val(meta), val(meta), tsv ]
    aracne_subnets = ARACNE3.out.subnets                // channel: [ val(meta), dir ]
    viper_results  = VIPER.out.viper_results            // channel: [ val(meta), tsv ]
    versions       = ch_versions                        // channel: [ versions.yml ]
}

