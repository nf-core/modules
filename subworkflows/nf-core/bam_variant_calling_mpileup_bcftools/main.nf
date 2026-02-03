include { BCFTOOLS_MPILEUP                 } from '../../../modules/nf-core/bcftools/mpileup'
include { BCFTOOLS_MERGE                   } from '../../../modules/nf-core/bcftools/merge'
include { BCFTOOLS_ANNOTATE                } from '../../../modules/nf-core/bcftools/annotate'

workflow BAM_VARIANT_CALLING_MPILEUP_BCFTOOLS {

    take:
    ch_bam     // channel: [ [id], bam, bai ]
    ch_posfile // channel: [ [panel_id, chr], posfile_comma]
    ch_fasta   // channel: [ [genome], fasta, fai]
    merge_key  // val    : [ "id" ]

    main:
    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

    ch_mpileup = ch_bam
        .combine(ch_posfile)
        .map{metaI, bam, _bai, metaPC, tsv ->
                [metaI + metaPC, bam, tsv]
        }

    BCFTOOLS_MPILEUP(
        ch_mpileup,
        ch_fasta,
        false
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_MPILEUP.out.stats.map{ it -> it[1] })

    // Branch depending on number of files
    ch_all_vcf = BCFTOOLS_MPILEUP.out.vcf
        .join(BCFTOOLS_MPILEUP.out.tbi)
        .map{ metaIPC, vcf, tbi -> // Get all keys except merge_key
            def groupKeys = metaIPC.keySet().findAll { it != merge_key }
            def groupMeta = metaIPC.subMap(groupKeys)
            [groupMeta, [metaIPC, vcf, tbi]]
        }
        .groupTuple(sort: { it1, it2 -> it1[0][merge_key] <=> it2[0][merge_key] }) // Sort by id
        .map{ metaPC, filestups ->
            // Create new meta with merge_key set to merged_value
            def newMeta = metaPC + [(merge_key): merged_value, metas: filestups.collect{it -> it[0]}]
            [
                newMeta,
                filestups.collect{it -> it[1]},
                filestups.collect{it -> it[2]},
                filestups.collect{it -> it[1]}.size()
            ]
        } // Compute number of records
        .branch{it ->
            one: it[3] == 1
            more: it[3] > 1
        }

    // Merge VCFs all individuals
    BCFTOOLS_MERGE(
        ch_all_vcf.more.map{it -> [it[0], it[1], it[2], []] },
        ch_fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions.first())

    // Mix all vcfs
    ch_to_concat = ch_all_vcf.one
        .map{it -> [it[0]["metas"][0], it[1][0], it[2][0]] }
        .mix(
            BCFTOOLS_MERGE.out.vcf
                .join(BCFTOOLS_MERGE.out.tbi.mix(
                    BCFTOOLS_MERGE.out.csi
                ))
        )

    // Mix all vcfs
    ch_to_concat = ch_all_vcf.one
        .map{it -> [it[0]["metas"][0], it[1][0], it[2][0]] }
        .mix(
            BCFTOOLS_MERGE.out.vcf
                .join(BCFTOOLS_MERGE.out.tbi.mix(
                    BCFTOOLS_MERGE.out.csi
                ))
        )

    // Merge all chromosomes
    VCF_CONCATENATE_BCFTOOLS(ch_to_concat)
    ch_versions = ch_versions.mix(VCF_CONCATENATE_BCFTOOLS.out.versions.first())

    // Annotate the variants
    BCFTOOLS_ANNOTATE(VCF_CONCATENATE_BCFTOOLS.out.vcf_index
        .combine(channel.of([[], [], [], []]))
    )
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())

    // Output
    ch_output = BCFTOOLS_ANNOTATE.out.vcf
        .join(BCFTOOLS_ANNOTATE.out.tbi.mix(
            BCFTOOLS_ANNOTATE.out.csi
        ))
        .map{ metaIPC, vcf, index -> [metaIPC + [ variantcaller:'bcftools' ], vcf, index] }


    emit:
    vcf_index     = ch_output        // channel: [ [id, panel], vcf, index ]
    versions      = ch_versions      // channel: [ versions.yml ]
    multiqc_files = ch_multiqc_files
}
