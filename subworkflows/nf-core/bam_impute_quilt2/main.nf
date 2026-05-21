include { QUILT_QUILT2    } from '../../../modules/nf-core/quilt/quilt2/main'
include { GLIMPSE2_LIGATE } from '../../../modules/nf-core/glimpse2/ligate/main'
include { BCFTOOLS_INDEX  } from '../../../modules/nf-core/bcftools/index/main'

workflow BAM_IMPUTE_QUILT2 {
    take:
    ch_input           // channel (mandatory): [ meta, bam/cram, bai/crai, bampaths, bamnames ]
    ch_reference_panel // channel (mandatory): [ meta, reference_vcf, reference_index ]
    ch_chunks          // channel (optional):  [ meta, chr, start, end ]
    ch_map             // channel (optional):  [ meta, genetic_map ]
    ch_fasta           // channel (optional):  [ meta, fasta, fai ]
    n_gen              // integer: Number of generations since founding or mixing
    buffer             // integer: Buffer of region to perform imputation over

    main:

    // Make final channel with parameters
    ch_chunks_counts = ch_chunks
        .groupTuple()
        .map { metaPC, chr, _start, _end ->
            [metaPC, chr.size()]
        }

    ch_parameters = ch_reference_panel
        .combine(ch_map, by: 0)
        .combine(ch_chunks, by: 0)
        .combine(ch_chunks_counts, by: 0)

    ch_parameters.ifEmpty {
        error("ERROR: join operation resulted in an empty channel. Please provide a valid ch_chunks and ch_map channel as input.")
    }

    ch_bam_params = ch_input
        .combine(ch_parameters)
        .map { meta_input, bam, bai, bampath, bamname, meta_panel, reference_vcf, reference_index, genetic_map, chr, start, end, region_size ->
            def regionout = "${chr}"
            def regionoutPadded = "${chr}"
            if (start != [] && end != []) {
                def paddedStart = String.format('%010d', start as long)
                def paddedEnd = String.format('%010d', end as long)
                regionoutPadded = "${chr}:${paddedStart}-${paddedEnd}"
                regionout = "${chr}:${start}-${end}"
            }

            [
                meta_panel + meta_input + ["regionout": regionout, "regionoutPadded": regionoutPadded, "regionSize": region_size],
                bam, bai,
                bampath, bamname,
                reference_vcf, reference_index,
                [], [], [],
                chr, start, end,
                n_gen, buffer,
                genetic_map,
            ]
        }

    QUILT_QUILT2(ch_bam_params, ch_fasta)

    ligate_input = QUILT_QUILT2.out.vcf
        .join(QUILT_QUILT2.out.tbi)
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

    ch_vcf_index = GLIMPSE2_LIGATE.out.merged_variants.join(
        BCFTOOLS_INDEX.out.tbi.mix(BCFTOOLS_INDEX.out.csi),
        failOnMismatch: true,
        failOnDuplicate: true,
    )

    emit:
    vcf_index = ch_vcf_index // channel: [ meta, vcf, tbi/csi ]
}
