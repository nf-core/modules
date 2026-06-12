include { STITCH                                  } from '../../../modules/nf-core/stitch/main'
include { GLIMPSE2_LIGATE                         } from '../../../modules/nf-core/glimpse2/ligate/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_PHASE  } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_LIGATE } from '../../../modules/nf-core/bcftools/index/main'

workflow BAM_IMPUTE_STITCH {
    take:
    ch_input // channel (mandatory):   [ [id], [bam], [bai], bampaths, bamnames ]
    ch_posfile // channel (mandatory):   [ [panel, chr], posfile ]
    ch_chunks // channel (optional) :   [ [panel, chr], chr, start, end ]
    ch_map // channel (optional) :   [ [panel, chr], map ]
    ch_fasta // channel (optional) :   [ [genome], fa, fai ]
    k_val // integer:   k value for STITCH
    n_gen // integer:   number of generations for STITCH
    buffer // integer:   Buffer size around each chunk for STITCH
    seed // value  :   seed for random number generator

    main:

    // Make final channel with parameters
    ch_chunks_counts = ch_chunks
        .groupTuple()
        .map { metaPC, chr, _start, _end ->
            [metaPC, chr.size()]
        }

    ch_parameters = ch_posfile
        .combine(ch_map, by: 0)
        .combine(ch_chunks, by: 0)
        .combine(ch_chunks_counts, by: 0)

    ch_parameters.ifEmpty {
        error("ERROR: join operation resulted in an empty channel. Please provide a valid ch_chunks and ch_map channel as input.")
    }

    ch_bam_params = ch_input
        .combine(ch_parameters)
        .map { metaI, bam, bai, bampath, bamname, metaPC, posfile, gmap, chr, start, end, region_size ->
            if (!chr) {
                error("ERROR: chromosome is not provided in ch_chunks.")
            }
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
                posfile,
                [],
                gmap,
                [],
                chr, start, end, buffer,
                k_val, n_gen,
            ]
        }

    STITCH(ch_bam_params, ch_fasta, seed)

    // Index imputed annotated VCF
    BCFTOOLS_INDEX_PHASE(STITCH.out.vcf)

    // Ligate all phased files in one and index it
    ligate_input = STITCH.out.vcf
        .join(
            BCFTOOLS_INDEX_PHASE.out.index,
            failOnMismatch: true,
            failOnDuplicate: true,
        )
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

    BCFTOOLS_INDEX_LIGATE(GLIMPSE2_LIGATE.out.merged_variants)

    // Join imputed and index files
    ch_vcf_index = GLIMPSE2_LIGATE.out.merged_variants.join(
        BCFTOOLS_INDEX_LIGATE.out.index,
        failOnMismatch: true,
        failOnDuplicate: true,
    )

    emit:
    vcf_index = ch_vcf_index // channel:   [ [id, chr], vcf, tbi ]
}
