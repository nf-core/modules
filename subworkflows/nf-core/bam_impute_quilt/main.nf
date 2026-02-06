include { QUILT_QUILT     } from '../../../modules/nf-core/quilt/quilt'
include { GLIMPSE2_LIGATE } from '../../../modules/nf-core/glimpse2/ligate'
include { BCFTOOLS_INDEX  } from '../../../modules/nf-core/bcftools/index'

workflow BAM_IMPUTE_QUILT {
    take:
    ch_input // channel (mandatory):   [ [id], [bam], [bai], bampaths, bamnames ]
    ch_hap_legend // channel (mandatory):   [ [panel, chr], hap, legend ]
    ch_posfile // channel (mandatory):   [ [panel, chr], posfile ]
    ch_chunks // channel (optional) :   [ [panel, chr], chr, start, end ]
    ch_map // channel (optional) :   [ [panel, chr], map ]
    ch_fasta // channel (optional) :   [ [genome], fa, fai ]
    n_gen // integer: Number of generations since founding or mixing
    buffer // integer: Buffer of region to perform imputation over

    main:

    ch_versions = channel.empty()

    // Make final channel with parameters
    ch_parameters = ch_hap_legend
        .combine(ch_posfile, by: 0)
        .combine(ch_map, by: 0)
        .combine(ch_chunks, by: 0)

    ch_parameters.ifEmpty {
        error("ERROR: join operation resulted in an empty channel. Please provide a valid ch_chunks and ch_map channel as input.")
    }

    ch_bam_params = ch_input
        .combine(ch_parameters)
        .map { metaI, bam, bai, bampath, bamname, metaPC, hap, legend, posfile, gmap, chr, start, end ->
            def regionout = "${chr}"
            if (start != [] && end != []) {
                regionout = "${chr}:${start}-${end}"
            }
            [
                metaPC + metaI + ["regionout": regionout],
                bam,
                bai,
                bampath,
                bamname,
                hap,
                legend,
                posfile,
                [],
                [],
                chr,
                start,
                end,
                n_gen,
                buffer,
                gmap,
            ]
        }

    QUILT_QUILT(ch_bam_params, ch_fasta)
    ch_versions = ch_versions.mix(QUILT_QUILT.out.versions.first())

    // Ligate all phased files in one and index it
    ligate_input = QUILT_QUILT.out.vcf
        .join(QUILT_QUILT.out.tbi)
        .map { meta, vcf, index ->
            def keysToKeep = meta.keySet() - ['regionout']
            [meta.subMap(keysToKeep), vcf, index]
        }
        .groupTuple()

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
    versions  = ch_versions // channel:   [ versions.yml ]
}
