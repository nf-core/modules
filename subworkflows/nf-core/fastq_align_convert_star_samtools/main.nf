include { SAMTOOLS_CONVERT  } from '../../../modules/nf-core/samtools/convert/main'
include { STAR_ALIGN        } from '../../../modules/nf-core/star/align/main'
include { SAMTOOLS_INDEX    } from '../../../modules/nf-core/samtools/index/main.nf'

workflow FASTQ_ALIGN_CONVERT_STAR_SAMTOOLS {
    take:
    reads           // channel: [ meta, [ fastq1, fastq2 ]]
    index           // channel: [ meta, index ]
    gtf             // channel: [ meta, gtf ]
    fasta           // channel: [ meta, fasta ]
    fai             // channel: [ meta, fai ]
    ignore_sjdbgtf  // boolean
    cram            // boolean: Create CRAM files

    main:
    def ch_versions = Channel.empty()
    STAR_ALIGN(
        reads,
        index,
        gtf,
        ignore_sjdbgtf,
        "", // seq_platform is handled in the config
        ""  // seq_center is handled in the config
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    SAMTOOLS_INDEX(
        STAR_ALIGN.out.bam_sorted_aligned
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_bam_bai = STAR_ALIGN.out.bam_sorted_aligned
        .join(SAMTOOLS_INDEX.out.bai, failOnMismatch:true, failOnDuplicate:true)

    def ch_cram_crai = Channel.empty()
    if(cram) {
        SAMTOOLS_CONVERT(
            STAR_ALIGN.out.bam_sorted_aligned.map { meta, bam -> [ meta, bam, []] },
            fasta,
            fai
        )
        ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions.first())
        ch_cram_crai = SAMTOOLS_CONVERT.out.cram
            .join(SAMTOOLS_CONVERT.out.crai, failOnMismatch:true, failOnDuplicate:true)
    }

    emit:
    versions      = ch_versions                      // channel: [ versions ]
    bam_bai       = ch_bam_bai                       // channel: [ val(meta), path(bam), path(bai) ]
    cram_crai     = ch_cram_crai                     // channel: [ val(meta), path(cram), path(crai) ]
    junctions     = STAR_ALIGN.out.junction          // channel: [ val(meta), path(junction) ]
    spl_junc_tabs = STAR_ALIGN.out.spl_junc_tab      // channel: [ val(meta), path(spl_jun_tab) ]
    log_final     = STAR_ALIGN.out.log_final         // channel: [ val(meta), path(log_final) ]
    gene_count    = STAR_ALIGN.out.read_per_gene_tab // channel: [ val(meta), path(read_per_gene_tab) ]
}
