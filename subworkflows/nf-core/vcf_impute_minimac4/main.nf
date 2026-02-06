include { MINIMAC4_COMPRESSREF                    } from '../../../modules/nf-core/minimac4/compressref'
include { MINIMAC4_IMPUTE                         } from '../../../modules/nf-core/minimac4/impute'
include { GLIMPSE2_LIGATE                         } from '../../../modules/nf-core/glimpse2/ligate'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_PHASE  } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_LIGATE } from '../../../modules/nf-core/bcftools/index'

workflow VCF_IMPUTE_MINIMAC4 {
    take:
    ch_input   // channel (mandatory): [ [id, chr], vcf, tbi ]
    ch_panel   // channel (mandatory): [ [panel, chr], vcf, tbi ]
    ch_posfile // channel (optional) : [ [panel, chr], sites_vcf, sites_index ]
    ch_chunks  // channel (optional) : [ [panel, chr], regionout ]
    ch_map     // channel (optional) : [ [panel, chr], map]

    main:

    ch_panel_branched = ch_panel.branch { _meta, file, _index ->
        def name = file.toString()
        vcf: name.matches(/.*\.(vcf|bcf)(\.gz)?$/)
        msav: name.endsWith('.msav')
        other: true
    }

    ch_panel_branched.other.map {
        error("ERROR: ch_panel files must be either VCF/BCF or MSAV.")
    }

    // Compress reference panel to MSAV format
    MINIMAC4_COMPRESSREF(ch_panel_branched.vcf)

    ch_panel_msav = MINIMAC4_COMPRESSREF.out.msav.mix(
        ch_panel_branched.msav.map { meta, file, _index -> [meta, file] }
    )

    // Channel with all reference and chunks informations
    // [ meta, reference panel msav, sites to impute vcf, sites index, region to impute, genetic map ]
    ch_panel_impute = ch_panel_msav
        .combine(ch_posfile, by: 0)
        .combine(ch_chunks, by: 0)
        .combine(ch_map, by: 0)

    ch_panel_impute.ifEmpty {
        error("ERROR: join operation resulted in an empty channel. Please provide a valid ch_posfile, ch_chunks and ch_map channel as input.")
    }

    // Prepare input channels for MINIMAC4
    ch_minimac4_input = ch_input
        .combine(ch_panel_impute)
        .map { metaI, target_vcf, target_tbi, metaPC, ref_msav, sites_vcf, sites_index, regionout, map ->
            [
                metaPC + metaI + ["regionout": regionout],
                target_vcf, target_tbi,
                ref_msav,
                sites_vcf, sites_index,
                map,
                regionout,
            ]
        }
    // Perform imputation
    MINIMAC4_IMPUTE(ch_minimac4_input)

    // Index the output VCF file
    BCFTOOLS_INDEX_PHASE(MINIMAC4_IMPUTE.out.vcf)

    // Ligate all phased files in one and index it
    ligate_input = MINIMAC4_IMPUTE.out.vcf
        .join(
            BCFTOOLS_INDEX_PHASE.out.tbi.mix(BCFTOOLS_INDEX_PHASE.out.csi),
            failOnMismatch: true,
            failOnDuplicate: true,
        )
        .map { meta, vcf, index ->
            def keysToKeep = meta.keySet() - ['regionout']
            [meta.subMap(keysToKeep), vcf, index]
        }
        .groupTuple()

    GLIMPSE2_LIGATE(ligate_input)

    BCFTOOLS_INDEX_LIGATE(GLIMPSE2_LIGATE.out.merged_variants)

    // Join imputed and index files
    ch_vcf_index = GLIMPSE2_LIGATE.out.merged_variants.join(
        BCFTOOLS_INDEX_LIGATE.out.tbi.mix(BCFTOOLS_INDEX_LIGATE.out.csi),
        failOnMismatch: true,
        failOnDuplicate: true,
    )

    emit:
    vcf_index = ch_vcf_index // channel: [ [id, panel, chr], vcf, index ]
}
