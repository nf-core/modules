include { CONCOCT_CUTUPFASTA           } from '../../../modules/nf-core/concoct/cutupfasta/main.nf'
include { CONCOCT_CONCOCTCOVERAGETABLE } from '../../../modules/nf-core/concoct/concoctcoveragetable/main.nf'
include { CONCOCT_CONCOCT              } from '../../../modules/nf-core/concoct/concoct/main.nf'
include { CONCOCT_MERGECUTUPCLUSTERING } from '../../../modules/nf-core/concoct/mergecutupclustering/main.nf'
include { CONCOCT_EXTRACTFASTABINS     } from '../../../modules/nf-core/concoct/extractfastabins/main.nf'

workflow FASTA_BINNING_CONCOCT {

    take:
    ch_fasta // channel (mandatory): [ val(meta), [ fasta ] ] (raw contigs from assembly)
    ch_bam   // channel (mandatory): [ val(meta), [ bam ], [bai]] (bam files of original FASTQ Files mapped back to each contig. meta must correspond to ch_fasta)

    main:
    ch_versions = Channel.empty()

    // required to create bedfile due to coverage table
    produce_bedfile = true

    CONCOCT_CUTUPFASTA ( ch_fasta, produce_bedfile )
    ch_versions = ch_versions.mix(CONCOCT_CUTUPFASTA.out.versions.first())

    ch_cutupfasta_for_concoctcoveragetable = CONCOCT_CUTUPFASTA.out.bed
                                                .join( ch_bam, failOnMismatch: true )

    CONCOCT_CONCOCTCOVERAGETABLE ( ch_cutupfasta_for_concoctcoveragetable )
    ch_versions = ch_versions.mix(CONCOCT_CONCOCTCOVERAGETABLE.out.versions.first())

    ch_concoctcoveragetable_for_concoctconcoct = CONCOCT_CONCOCTCOVERAGETABLE.out.tsv
                                                    .join(CONCOCT_CUTUPFASTA.out.fasta, failOnMismatch: true)

    CONCOCT_CONCOCT( ch_concoctcoveragetable_for_concoctconcoct )
    ch_versions = ch_versions.mix(CONCOCT_CONCOCT.out.versions.first())

    CONCOCT_MERGECUTUPCLUSTERING ( CONCOCT_CONCOCT.out.clustering_csv )
    ch_versions = ch_versions.mix( CONCOCT_MERGECUTUPCLUSTERING.out.versions.first())

    ch_mergecutupclustering_for_extractfastabins = ch_fasta
                                                    .join(CONCOCT_MERGECUTUPCLUSTERING.out.csv, failOnMismatch: false)

    CONCOCT_EXTRACTFASTABINS ( ch_mergecutupclustering_for_extractfastabins )
    ch_versions = ch_versions.mix(CONCOCT_EXTRACTFASTABINS.out.versions.first())

    emit:
    coverage_table      = CONCOCT_CONCOCTCOVERAGETABLE.out.tsv     // channel: [ val(meta), [ tsv ] ]

    original_csv        = CONCOCT_CONCOCT.out.original_data_csv    // channel: [ val(meta), [ csv ] ]
    raw_clustering_csv  = CONCOCT_CONCOCT.out.clustering_csv       // channel: [ val(meta), [ csv ] ]
    pca_original        = CONCOCT_CONCOCT.out.pca_components_csv   // channel: [ val(meta), [ csv ] ]
    pca_transformed     = CONCOCT_CONCOCT.out.pca_transformed_csv  // channel: [ val(meta), [ csv ] ]

    cluster_table       = CONCOCT_MERGECUTUPCLUSTERING.out.csv     // channel: [ val(meta), [ csv ] ]
    bins                = CONCOCT_EXTRACTFASTABINS.out.fasta       // channel: [ val(meta), [ fasta ] ]

    versions = ch_versions                                         // channel: [ versions.yml ]
}
