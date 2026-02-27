include { BCFTOOLS_MPILEUP                 } from '../../../modules/nf-core/bcftools/mpileup'
include { BCFTOOLS_MERGE                   } from '../../../modules/nf-core/bcftools/merge'
include { BCFTOOLS_ANNOTATE                } from '../../../modules/nf-core/bcftools/annotate'
include { VCF_GATHER_BCFTOOLS              } from '../../../subworkflows/nf-core/vcf_gather_bcftools'

workflow BAM_VARIANT_CALLING_MPILEUP_BCFTOOLS {

    take:
    ch_bam                   // channel: [ [id], bam, bai ]
    ch_posfile               // channel: [ [panel_id, chr], posfile_comma]
    ch_fasta                 // channel: [ [genome], fasta, fai ]
    meta_sample_merge_key    // val    : [ "id" ]
    meta_sample_merge_value  // val    : [ "all_samples" ]
    meta_region_gather_keys  // val    : [ "panel_id", "id" ]
    sort_region_gather       // val    : boolean
    annotate                 // val    : boolean

    main:
    ch_multiqc_files = channel.empty()

    ch_mpileup = ch_bam
        .combine(ch_posfile)
        .map{metaI, bam, _bai, metaPC, tsv ->
                [metaI + metaPC, bam, tsv, tsv]
        }

    def posfile_count = ch_posfile
        .map{ _meta, posfile -> posfile}
        .collect()
        .map { posfile -> posfile.size() }

    BCFTOOLS_MPILEUP(
        ch_mpileup,
        ch_fasta,
        false
    )
    ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_MPILEUP.out.stats.map{ _meta, stats -> stats })

    // Branch depending on number of files
    ch_all_vcf = BCFTOOLS_MPILEUP.out.vcf
        .join(BCFTOOLS_MPILEUP.out.tbi)
        .map{ metaIPC, vcf, tbi -> // Get all keys except merge_key
            def groupKeys = metaIPC.keySet().findAll { key -> key != meta_sample_merge_key }
            def groupMeta = metaIPC.subMap(groupKeys)
            [groupMeta, [metaIPC, vcf, tbi]]
        }
        .groupTuple(sort: { it1, it2 -> it1[0][meta_sample_merge_key] <=> it2[0][meta_sample_merge_key] }) // Sort by id
        .map{ metaPC, filestups ->
            // Create new meta with meta_sample_merge_key set to meta_sample_merge_value
            def newMeta = metaPC + [
                (meta_sample_merge_key): meta_sample_merge_value,
                metas: filestups.collect{ meta, _vcf, _index -> meta }
            ]
            [
                newMeta,
                filestups.collect{_meta, vcf, _index -> vcf},
                filestups.collect{_meta, _vcf, index -> index},
                filestups.collect{_meta, vcf, _index -> vcf}.size()
            ]
        }
        .branch{ _meta, _vcf, _index, size ->
            one: size == 1
            more: size > 1
        }

    // Merge VCFs all individuals
    BCFTOOLS_MERGE(
        ch_all_vcf.more.map{
            meta, vcf_list, index_list, _size -> [ meta, vcf_list, index_list, [] ]
        },
        ch_fasta
    )

    // Mix all vcfs
    ch_to_concat = ch_all_vcf.one
        .map{ meta, vcf_list, index_list, _size -> [
            meta, vcf_list[0], index_list[0]
        ] }
        .mix(
            BCFTOOLS_MERGE.out.vcf
                .join(BCFTOOLS_MERGE.out.index)
        )

    // Merge all chromosomes
    VCF_GATHER_BCFTOOLS(
        ch_to_concat.combine(posfile_count),
        meta_region_gather_keys,
        sort_region_gather
    )

    if (annotate) {
        // Annotate the variants
        BCFTOOLS_ANNOTATE(VCF_GATHER_BCFTOOLS.out.vcf_index
            .combine(channel.of([[], [], [], [], []]))
        )
        // Output
        ch_output = BCFTOOLS_ANNOTATE.out.vcf
            .join(BCFTOOLS_ANNOTATE.out.tbi.mix(
                BCFTOOLS_ANNOTATE.out.csi
            ))
    } else {
        // Output without annotation
        ch_output = VCF_GATHER_BCFTOOLS.out.vcf_index
    }

    emit:
    vcf_index     = ch_output        // channel: [ [id, panel], vcf, index ]
    multiqc_files = ch_multiqc_files
}
