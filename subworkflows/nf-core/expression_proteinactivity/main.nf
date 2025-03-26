//
// Infer protein activity interactions of a set of
//     gene regulators from gene expression data
//

include { SALMON_INDEX } from '../../../modules/nf-core/salmon/index/main'
include { SALMON_QUANT } from '../../../modules/nf-core/salmon/quant/main'
// include { ARACNE3 } from '../../../modules/nf-core/aracne3/main'
// include { VIPER } from '../../../modules/nf-core/viper/main'

workflow EXPRESSION_PROTEINACTIVITY {

    take:
    // ch_reads             // channel: [ val(meta), [ reads ] ]
    // ch_genome_fasta      // channel: /path/to/genome.fasta
    // ch_transcript_fasta  // channel: /path/to/transcript.fasta
    // ch_gtf               // channel: /path/to/genome.gtf
    // ch_index             // channel: /path/to/salmon/index/
    // make_index           // boolean: Whether to create salmon index before running salmon quant
    ch_regulators        // channel: /path/to/regulators_genes.txt
    ch_expression_matric // channel: /path/to/expression_quant.sf
    ch_network           // channel: /path/to/regulon_network.txt
    ch_sample_expression // channel: [ val(meta), [ reads ] ]

    main:

    ch_versions = Channel.empty()

    //
    // Generate Gene Expression Quantification
    //

    // // Create Salmon index if required
    // if (make_index) {
    //     ch_index = SALMON_INDEX ( ch_genome_fasta, ch_transcript_fasta ).index
    //     ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
    // }

    // // Salmon Quantification
    // def lib_type = 'A'
    // def alignment_mode = false
    // SALMON_QUANT ( ch_reads, ch_index, ch_gtf, ch_transcript_fasta, alignment_mode, lib_type )
    // ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())

    // Run ARACNE
    ARACNE3 ( ch_expression_, ch_regulators ) // Get tsv from results_dir

    // Run VIPER
    VIPER ( ARACNE3.out.consensus_network, ch_sample_expression )

    emit:
    // index             = ch_index                           // channel: [ index ]
    // results           = SALMON_QUANT.out.results           // channel: [ val(meta), results_dir ]
    // json_info         = SALMON_QUANT.out.json_info         // channel: [ val(meta), json_info ]
    // lib_format_counts = SALMON_QUANT.out.lib_format_counts // channel: [ val(meta), lib_format_counts ]

    consensus_network = ARACNE3.out.consensus_network      // channel: [ tsv ]
    subnets           = ARACNE3.out.subnets                 // channel: [ dir ]

    consensus_network = ARACNE.out.consensus_network      // channel: [ tsv ]
    subnets           = ARACNE.out.subnets                 // channel: [ dir ]

    versions          = ch_versions                        // channel: [ versions.yml ]
}


