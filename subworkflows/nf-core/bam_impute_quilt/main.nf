include { QUILT_QUILT     } from '../../../modules/nf-core/quilt/quilt'
include { GLIMPSE2_LIGATE } from '../../../modules/nf-core/glimpse2/ligate'
include { BCFTOOLS_INDEX  } from '../../../modules/nf-core/bcftools/index'

workflow BAM_IMPUTE_QUILT {
    take:
    ch_input // channel (mandatory):   [ [id], [bam], [bai], bampaths, bamnames ]
    ch_hap_legend_posfile // channel (mandatory):   [ [panel, chr], hap, legend, posfile ]
    ch_chunks // channel (optional) :   [ [panel, chr], chr, start, end ]
    ch_map // channel (optional) :   [ [panel, chr], map ]
    ch_fasta // channel (optional) :   [ [genome], fa, fai ]
    n_gen // integer: Number of generations since founding or mixing
    buffer // integer: Buffer of region to perform imputation over

    main:

    // Make final channel with parameters
    ch_chunks_counts = ch_chunks
        .groupTuple()
        .map { metaPC, chr, _start, _end ->
            [metaPC, chr.size()]
        }

    ch_parameters = ch_hap_legend_posfile
        .combine(ch_map, by: 0)
        .combine(ch_chunks, by: 0)
        .combine(ch_chunks_counts, by: 0)

    ch_parameters.ifEmpty {
        error("ERROR: join operation resulted in an empty channel. Please provide a valid ch_chunks and ch_map channel as input.")
    }

    ch_bam_params = ch_input
        .combine(ch_parameters)
        .map { metaI, bam, bai, bampath, bamname, metaPC, hap, legend, posfile, gmap, chr, start, end, region_size ->
            def regionout = "${chr}"
            def regionoutPadded = "${chr}"
            if (start != [] && end != []) {
                def paddedStart = String.format('%010d', start as long)
                def paddedEnd = String.format('%010d', end as long)
                regionoutPadded = "${chr}:${paddedStart}-${paddedEnd}"
                regionout = "${chr}:${start}-${end}"
            }
            [
                metaPC + metaI + ["regionout": regionout, "regionoutPadded": regionoutPadded, "regionSize": region_size],
                bam, bai,
                bampath, bamname,
                hap, legend,
                posfile,
                [], [],
                chr, start, end,
                n_gen, buffer,
                gmap,
            ]
        }

    QUILT_QUILT(ch_bam_params, ch_fasta)

    // Ligate all phased files in one and index it
    ligate_input = QUILT_QUILT.out.vcf
        .join(QUILT_QUILT.out.tbi)
        .map { meta, vcf, index ->
            def keysToKeep = meta.keySet() - ['regionout', 'regionoutPadded', 'regionSize']
            [
                groupKey(meta.subMap(keysToKeep), meta.regionSize),
                vcf, index,
            ]
        }
        .groupTuple()
        .map { groupKeyObj, vcf, index ->
            // Extract the actual meta from the groupKey
            def meta = groupKeyObj.getGroupTarget()
            [meta, vcf, index]
        }

    GLIMPSE2_LIGATE(ligate_input)

    BCFTOOLS_INDEX(GLIMPSE2_LIGATE.out.merged_variants)

    // Join imputed and index files
    ch_vcf_index = GLIMPSE2_LIGATE.out.merged_variants.join(
        BCFTOOLS_INDEX.out.tbi.mix(BCFTOOLS_INDEX.out.csi),
        failOnMismatch: true,
        failOnDuplicate: true,
    )

    emit:
    vcf_index = ch_vcf_index // channel:   [ [id, chr], vcf, tbi ]
}
