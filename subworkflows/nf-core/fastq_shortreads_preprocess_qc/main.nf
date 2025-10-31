// statistics
include { FASTQ_GENERATE_STATISTICS } from '../../../subworkflows/nf-core/fastq_generate_statistics/main.nf'
// preprocess
include { FASTQ_PREPROCESS          } from '../../../subworkflows/nf-core/fastq_preprocess/main.nf' 
/*
// preprocessing
include { SEQKIT_SEQ            } from '../../../modules/nf-core/seqkit/seq/main'
// barcoding
include { UMITOOLS_EXTRACT      } from '../../../modules/nf-core/umitools/extract/main'
// adapter removal and merging
include { FASTP                 } from '../../../modules/nf-core/fastp/main'
// complexity filtering
include { PRINSEQPLUSPLUS       } from '../../../modules/nf-core/prinseqplusplus/main'
// deduplication
include { BBMAP_CLUMPIFY        } from '../../../modules/nf-core/bbmap/clumpify/main'
*/
// host decontamination
//include { HOSTILE_FETCH         } from '../../../modules/nf-core/hostile/fetch/main'
//include { HOSTILE_CLEAN         } from '../../../modules/nf-core/hostile/clean/main'// final concatenation
//include { CAT_FASTQ                         } from '../../../modules/nf-core/cat/fastq/main'

workflow FASTQ_SHORTREADS_PREPROCESS_QC {

    take:
    ch_reads          // channel: [ val(meta), [ fastq ] ]
    skip_fastqc       // boolean
    skip_seqfu_check  // boolean
    skip_seqfu_stats  // boolean
    skip_seqkit_stats // boolean
    skip_seqtk_comp   // boolean

    skip_seqkit_sana_pair // boolean
    skip_seqkit_seq       // boolean
    skip_seqkit_replace   // boolean
    skip_seqkit_rmdup     // boolean


    main:
    ch_versions = Channel.empty()

    FASTQ_GENERATE_STATISTICS(
        ch_reads,
        skip_fastqc,       
        skip_seqfu_check, 
        skip_seqfu_stats,  
        skip_seqkit_stats,
        skip_seqtk_comp
    )
    ch_versions    = ch_versions.mix(FASTQ_GENERATE_STATISTICS.out.versions)


    FASTQ_PREPROCESS(
        ch_reads,
        skip_seqkit_sana_pair,
        skip_seqkit_seq,
        skip_seqkit_replace,
        skip_seqkit_rmdup
    )

    ch_reads = FASTQ_PREPROCESS.out.reads
    ch_versions    = ch_versions.mix(FASTQ_PREPROCESS.out.versions)
/*
    // barcoding
    // umi_reads is defined as ch_reads in case this process is skipped
    umi_reads = ch_reads
    umi_log = Channel.empty()
    // if process is executed, umi_reads will take a new ouptut
    if (!skip_umitools_extract) {
        UMITOOLS_EXTRACT( ch_reads )
        umi_reads = UMITOOLS_EXTRACT.out.reads
        umi_log = UMITOOLS_EXTRACT.out.log
        ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())

        // Discard R1 / R2 if required
        if (umi_discard_read in [1, 2]) {
            // if data is single end, remove paired reads, but how are these generated?
            UMITOOLS_EXTRACT.out.reads
                .map { meta, reads ->
                    meta.single_end ? [meta, reads] : [meta + ['single_end': true], reads[umi_discard_read % 2]]
                }
                .set { umi_reads }
        }
    }

    // adapter removal and merging
    ch_reads_preprocessed = umi_reads
    if (!skip_fastp){
        // FASTP requirements
            //tuple val(meta), path(reads), path(adapter_fasta)
            //val   discard_trimmed_pass
            //val   save_trimmed_fail
            //val   save_merged
        
        // EXAMPLE
            //                 adapter_fasta        = [] // empty list for no adapter file!
            //    discard_trimmed_pass = false
            //    save_trimmed_fail    = false
            //    save_merged          = false

        fastp_input = umi_reads.combine( fastp_adapter_fasta )
        FASTP(
            fastp_input,
            fastp_discard_trimmed_pass,
            fastp_save_trimmed_fail,
            fastp_save_merged )

        FASTP.out.reads.set { ch_reads_preprocessed }

        ch_versions = ch_versions.mix(FASTP.out.versions.first())

    }

    // complexity filtering
    ch_reads_complexity_filtered = ch_reads_preprocessed
    if (!skip_complexity_filtering) {
        PRINSEQPLUSPLUS( ch_reads_preprocessed )
        ch_versions = ch_versions.mix(PRINSEQPLUSPLUS.out.versions.first())
        ch_reads_complexity_filtered = PRINSEQPLUSPLUS.out.good_reads
    }

    // deduplication
    ch_reads_deduplicated = ch_reads_complexity_filtered
    if (!skip_deduplication) {
        BBMAP_CLUMPIFY( ch_reads_complexity_filtered )
        ch_versions = ch_versions.mix(BBMAP_CLUMPIFY.out.versions.first())
        ch_reads_deduplicated = BBMAP_CLUMPIFY.out.reads
    }

    // host decontamination

    // final concatenation
    // TODO
    // if (!skip_final_concatenation) {
    //     CAT_FASTQ( ... )
    //     ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())
    // }

    // post - statistics
    // TODO
    // ...

*/

    emit:
    fastqc_html   = FASTQ_GENERATE_STATISTICS.out.fastqc_html
    fastqc_zip    = FASTQ_GENERATE_STATISTICS.out.fastqc_zip
    seqfu_check   = FASTQ_GENERATE_STATISTICS.out.seqfu_check
    seqfu_stats   = FASTQ_GENERATE_STATISTICS.out.seqfu_stats
    seqfu_multiqc = FASTQ_GENERATE_STATISTICS.out.seqfu_multiqc
    seqkit_stats  = FASTQ_GENERATE_STATISTICS.out.seqkit_stats
    seqtk_stats   = FASTQ_GENERATE_STATISTICS.out.seqtk_stats
    versions      = ch_versions

}