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

    ch_versions = channel.empty()

    ch_parameters = ch_reference_panel
        .combine(ch_map, by: 0)
        .combine(ch_chunks, by: 0)

    ch_parameters.ifEmpty {
        error("ERROR: join operation resulted in an empty channel. Please provide a valid ch_chunks and ch_map channel as input.")
    }

    ch_bam_params = ch_input
        .combine(ch_parameters)
        .map { meta_input, bam, bai, bampath, bamname, meta_panel, reference_vcf, reference_index, genetic_map, chr, start, end ->
            def regionout = "${chr}"
            if (start != [] && end != []) {
                regionout = "${chr}:${start}-${end}"
            }

            [
                meta_panel + meta_input + [regionout: regionout],
                bam,
                bai,
                bampath,
                bamname,
                reference_vcf,
                reference_index,
                [],
                [],
                [],
                chr,
                start,
                end,
                n_gen,
                buffer,
                genetic_map,
            ]
        }

    QUILT_QUILT2(ch_bam_params, ch_fasta)
    ch_versions = ch_versions.mix(QUILT_QUILT2.out.versions_r_quilt)
    ch_versions = ch_versions.mix(QUILT_QUILT2.out.versions_r_base)

    ligate_input = QUILT_QUILT2.out.vcf
        .join(QUILT_QUILT2.out.tbi)
        .map { meta, vcf, index ->
            def keys_to_keep = meta.keySet() - ['regionout']
            [meta.subMap(keys_to_keep), vcf, index]
        }
        .groupTuple()

    GLIMPSE2_LIGATE(ligate_input)

    BCFTOOLS_INDEX(GLIMPSE2_LIGATE.out.merged_variants)

    ch_vcf_index = GLIMPSE2_LIGATE.out.merged_variants.join(
        BCFTOOLS_INDEX.out.tbi.mix(BCFTOOLS_INDEX.out.csi),
        failOnMismatch: true,
        failOnDuplicate: true,
    )

    emit:
    vcf_index = ch_vcf_index // channel: [ meta, vcf, tbi/csi ]
    versions  = ch_versions  // channel: [ versions.yml ] or topic-based version tuples
}
