include { BEAGLE5_BEAGLE                     } from '../../../modules/nf-core/beagle5/beagle'
include { BCFTOOLS_VIEW                      } from '../../../modules/nf-core/bcftools/view'
include { GLIMPSE2_LIGATE                    } from '../../../modules/nf-core/glimpse2/ligate'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1 } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2 } from '../../../modules/nf-core/bcftools/index'



workflow VCF_IMPUTE_BEAGLE5 {

    take:
    ch_input  // channel (mandatory): [ [id], vcf, tbi ]
    ch_panel  // channel (mandatory): [ [panel, chr], vcf, tbi ]
    ch_chunks // channel (optional) : [ [panel, chr], regionout ]
    ch_map    // channel (optional) : [ [chr], map]

    main:
    
    ch_versions = channel.empty()

    // Branch input files based on format
    ch_input
        .branch { _meta, vcf, _tbi ->
            bcf: vcf.toString().contains('.bcf')
            vcf: vcf.toString().contains('.vcf')
            other: true
        }
        .set { ch_input_branched }
    
    ch_input_branched.other.map{
        error "ERROR: ch_input files must be in VCF or BCF format."
    }

    // Convert BCF to VCF if necessary
    BCFTOOLS_VIEW(
        ch_input_branched.bcf,
        [], [], []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())

    // Combine VCF files
    ch_ready_vcf = ch_input_branched.vcf
        .mix(BCFTOOLS_VIEW.out.vcf
            .join(
                BCFTOOLS_VIEW.out.csi
                    .mix(BCFTOOLS_VIEW.out.tbi)
            )
        )

    // Prepare input channels for BEAGLE5 by combining VCF, panel, and map files
    ch_panel_map = ch_panel
        .combine(ch_map, by: 0)
        .combine(ch_chunks, by: 0)

    ch_panel_map.ifEmpty{
        error "ERROR: join operation resulted in an empty channel. Please provide a valid ch_panel and ch_map channel as input."
    }

    ch_beagle_input = ch_ready_vcf
        .combine(ch_panel_map)
        .map { metaI, input_vcf, input_index, metaPC, panel_vcf, panel_index, map, regionout -> [
            metaI + metaPC + ["regionout": regionout],
            input_vcf, input_index,
            panel_vcf, panel_index,
            map, [], [], regionout
        ]}

    // Run BEAGLE5 imputation
    BEAGLE5_BEAGLE(ch_beagle_input)
    ch_versions = ch_versions.mix(BEAGLE5_BEAGLE.out.versions.first())

    // Index the imputed VCF files
    BCFTOOLS_INDEX_1(BEAGLE5_BEAGLE.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_1.out.versions.first())

    // Ligate all phased files in one and index it
    ligate_input = BEAGLE5_BEAGLE.out.vcf
        .join(
            BCFTOOLS_INDEX_1.out.tbi
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
    vcf_index  = ch_vcf_index // channel: [ [id, chr, tools], vcf, index ]
    versions   = ch_versions  // channel: [ versions.yml ]
}
