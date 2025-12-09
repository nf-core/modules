include { STITCH                             } from '../../../modules/nf-core/stitch/main'
include { GLIMPSE2_LIGATE                    } from '../../../modules/nf-core/glimpse2/ligate/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1 } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2 } from '../../../modules/nf-core/bcftools/index/main'

workflow BAM_IMPUTE_STITCH {

    take:
    ch_input        // channel (mandatory):   [ [id], [bam], [bai], bampaths, bamnames ]
    ch_posfile      // channel (mandatory):   [ [panel, chr], posfile ]
    ch_chunks       // channel (optional) :   [ [panel, chr], chr, start, end ]
    ch_map          // channel (optional) :   [ [panel, chr], map ]
    ch_fasta        // channel (optional) :   [ [genome], fa, fai ]
    k_val           // integer:   k value for STITCH
    n_gen           // integer:   number of generations for STITCH
    seed            // value  :   seed for random number generator

    main:

    ch_versions      = channel.empty()

    // Make final channel with parameters
    ch_parameters = ch_posfile
        .combine(ch_map, by: 0)
        .combine(ch_chunks, by: 0)

    ch_parameters.ifEmpty{
        error "ERROR: join operation resulted in an empty channel. Please provide a valid ch_chunks and ch_map channel as input."
    }

    ch_bam_params = ch_input
        .combine(ch_parameters)
        .map{
            metaI, bam, bai, bampath, bamname, metaPC, posfile, gmap, chr, start, end ->
            if (!chr) {
                error "ERROR: chromosome is not provided in ch_chunks."
            }
            def regionout = "${chr}"
            if (start != [] && end != []) {
                regionout = "${chr}:${start}-${end}"
            }
            [
                metaPC + metaI + ["regionout": regionout],
                bam, bai, bampath, bamname, posfile, [], gmap, [], chr, start, end, k_val, n_gen
            ]
        }

    STITCH( ch_bam_params, ch_fasta, seed )
    ch_versions = ch_versions.mix(STITCH.out.versions.first())

    // Index imputed annotated VCF
    BCFTOOLS_INDEX_1(STITCH.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_1.out.versions.first())

    // Ligate all phased files in one and index it
    ligate_input = STITCH.out.vcf
        .join( BCFTOOLS_INDEX_1.out.tbi
            .mix(BCFTOOLS_INDEX_1.out.csi)
        )
        .map{ meta, vcf, index ->
            def keysToKeep = meta.keySet() - ['regionout']
            [ meta.subMap(keysToKeep), vcf, index ]
        }
        .groupTuple()

    GLIMPSE2_LIGATE( ligate_input )
    ch_versions = ch_versions.mix( GLIMPSE2_LIGATE.out.versions.first() )

    BCFTOOLS_INDEX_2( GLIMPSE2_LIGATE.out.merged_variants )
    ch_versions = ch_versions.mix( BCFTOOLS_INDEX_2.out.versions.first() )

    // Join imputed and index files
    ch_vcf_index = GLIMPSE2_LIGATE.out.merged_variants
        .join(
            BCFTOOLS_INDEX_2.out.tbi
                .mix(BCFTOOLS_INDEX_2.out.csi)
        )

    emit:
    vcf_index  = ch_vcf_index  // channel:   [ [id, chr], vcf, tbi ]
    versions   = ch_versions   // channel:   [ versions.yml ]

}
