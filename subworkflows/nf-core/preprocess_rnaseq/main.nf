import groovy.json.JsonSlurper

include { BBMAP_BBSPLIT                   } from '../../../modules/nf-core/bbmap/bbsplit'
include { CAT_FASTQ                       } from '../../../modules/nf-core/cat/fastq/main'
include { SORTMERNA                       } from '../../../modules/nf-core/sortmerna/main'
include { SORTMERNA as SORTMERNA_INDEX    } from '../../../modules/nf-core/sortmerna/main'

include { FASTQ_SUBSAMPLE_FQ_SALMON        } from '../fastq_subsample_fq_salmon'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../fastq_fastqc_umitools_trimgalore'
include { FASTQ_FASTQC_UMITOOLS_FASTP      } from '../fastq_fastqc_umitools_fastp'

def pass_trimmed_reads = [:]

public static String getSalmonInferredStrandedness(json_file) {
    def lib_type = new JsonSlurper().parseText(json_file.text).get('library_types')[0]
    def strandedness = 'reverse'
    if (lib_type) {
        if (lib_type in ['U', 'IU']) {
            strandedness = 'unstranded'
        } else if (lib_type in ['SF', 'ISF']) {
            strandedness = 'forward'
        } else if (lib_type in ['SR', 'ISR']) {
            strandedness = 'reverse'
        }
    }
    return strandedness
}

//
// Create MultiQC tsv custom content from a list of values
//
public static String multiqcTsvFromList(tsv_data, header) {
    def tsv_string = ""
    if (tsv_data.size() > 0) {
        tsv_string += "${header.join('\t')}\n"
        tsv_string += tsv_data.join('\n')
    }
    return tsv_string
}

def deprecation_message = """
WARNING: This subworkflow has been deprecated. Please use
nf-core/modules/subworkflows/fastq_qc_trim_filter_setstrandedness

Reason:
Subworkflow naming rules were introduced necessitating this move, see
https://nf-co.re/docs/guidelines/components/subworkflows#name-format-of-subworkflow-files
"""

workflow PREPROCESS_RNASEQ {

    take:
    ch_reads             // channel: [ val(meta), [ reads ] ]
    ch_fasta             // channel: /path/to/genome.fasta
    ch_transcript_fasta  // channel: /path/to/transcript.fasta
    ch_gtf               // channel: /path/to/genome.gtf
    ch_salmon_index      // channel: /path/to/salmon/index/ (optional)
    ch_sortmerna_index   // channel: /path/to/sortmerna/index/ (optional)
    ch_bbsplit_index     // channel: /path/to/bbsplit/index/ (optional)
    ch_ribo_db           // channel: /path/to/ Text file containing paths to fasta files (one per line) that will be used to create the database for SortMeRNA. (optional)
    skip_bbsplit         // boolean: Skip BBSplit for removal of non-reference genome reads.
    skip_fastqc          // boolean: true/false
    skip_trimming        // boolean: true/false
    skip_umi_extract     // boolean: true/false
    make_salmon_index    // boolean: Whether to create salmon index before running salmon quant
    make_sortmerna_index // boolean: Whether to create a sortmerna index before running sortmerna
    trimmer              // string (enum): 'fastp' or 'trimgalore'
    min_trimmed_reads    // integer: > 0
    save_trimmed         // boolean: true/false
    remove_ribo_rna      // boolean: true/false: whether to run sortmerna to remove rrnas
    with_umi             // boolean: true/false: Enable UMI-based read deduplication.
    umi_discard_read     // integer: 0, 1 or 2

    main:

    assert false: deprecation_message

    ch_versions        = Channel.empty()
    ch_filtered_reads  = Channel.empty()
    ch_trim_read_count = Channel.empty()
    ch_multiqc_files   = Channel.empty()

    ch_reads
        .branch {
            meta, fastqs ->
                single  : fastqs.size() == 1
                    return [ meta, fastqs.flatten() ]
                multiple: fastqs.size() > 1
                    return [ meta, fastqs.flatten() ]
        }
        .set { ch_fastq }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_filtered_reads }

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with TrimGalore!
    //
    if (trimmer == 'trimgalore') {
        FASTQ_FASTQC_UMITOOLS_TRIMGALORE (
            ch_filtered_reads,
            skip_fastqc,
            with_umi,
            skip_umi_extract,
            skip_trimming,
            umi_discard_read,
            min_trimmed_reads
        )
        ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
        ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count

        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)
        ch_multiqc_files = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
            .mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip)
            .mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log)
            .mix(ch_multiqc_files)
    }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with fastp
    //
    if (trimmer == 'fastp') {
        FASTQ_FASTQC_UMITOOLS_FASTP (
            ch_filtered_reads,
            skip_fastqc,
            with_umi,
            skip_umi_extract,
            umi_discard_read,
            skip_trimming,
            [],
            save_trimmed,
            save_trimmed,
            min_trimmed_reads
        )
        ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
        ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_read_count

        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)
        ch_multiqc_files = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip
            .mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip)
            .mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json.map{tuple(it[0], [it[1]])})
            .mix(ch_multiqc_files)
    }

    //
    // Get list of samples that failed trimming threshold for MultiQC report
    //

    ch_trim_read_count
        .map {
            meta, num_reads ->
                pass_trimmed_reads[meta.id] = true
                if (num_reads <= min_trimmed_reads.toFloat()) {
                    pass_trimmed_reads[meta.id] = false
                    return [ "$meta.id\t$num_reads" ]
                }
        }
        .collect()
        .map {
            tsv_data ->
                def header = ["Sample", "Reads after trimming"]
                multiqcTsvFromList(tsv_data, header)
        }
        .set { ch_fail_trimming_multiqc }

    ch_multiqc_files = ch_multiqc_files
        .mix(
            ch_fail_trimming_multiqc.collectFile(name: 'fail_trimmed_samples_mqc.tsv')
        )

    //
    // MODULE: Remove genome contaminant reads
    //
    if (!skip_bbsplit) {
        BBMAP_BBSPLIT (
            ch_filtered_reads,
            ch_bbsplit_index,
            [],
            [ [], [] ],
            false
        )

        BBMAP_BBSPLIT.out.primary_fastq
            .set { ch_filtered_reads }

        ch_versions = ch_versions.mix(BBMAP_BBSPLIT.out.versions.first())
    }

    //
    // MODULE: Remove ribosomal RNA reads
    //
    if (remove_ribo_rna) {
        ch_sortmerna_fastas = Channel.from(ch_ribo_db.readLines())
            .map { row -> file(row, checkIfExists: true) }
            .collect()
            .map{ ['rrna_refs', it] }

        if (make_sortmerna_index) {
            SORTMERNA_INDEX (
                [[],[]],
                ch_sortmerna_fastas,
                [[],[]]
            )
            ch_sortmerna_index = SORTMERNA_INDEX.out.index.first()
        }

        SORTMERNA (
            ch_filtered_reads,
            ch_sortmerna_fastas,
            ch_sortmerna_index
        )

        SORTMERNA.out.reads
            .set { ch_filtered_reads }

        ch_multiqc_files = ch_multiqc_files
            .mix(SORTMERNA.out.log)

        ch_versions = ch_versions.mix(SORTMERNA.out.versions.first())
    }

    // Branch FastQ channels if 'auto' specified to infer strandedness
    ch_filtered_reads
        .branch {
            meta, fastq ->
                auto_strand : meta.strandedness == 'auto'
                    return [ meta, fastq ]
                known_strand: meta.strandedness != 'auto'
                    return [ meta, fastq ]
        }
        .set { ch_strand_fastq }

    //
    // SUBWORKFLOW: Sub-sample FastQ files and pseudoalign with Salmon to auto-infer strandedness
    //
    // Return empty channel if ch_strand_fastq.auto_strand is empty so salmon index isn't created

    ch_fasta
        .combine(ch_strand_fastq.auto_strand)
        .map { it.first() }
        .first()
        .set { ch_genome_fasta }

    FASTQ_SUBSAMPLE_FQ_SALMON (
        ch_strand_fastq.auto_strand,
        ch_genome_fasta,
        ch_transcript_fasta,
        ch_gtf,
        ch_salmon_index,
        make_salmon_index
    )
    ch_versions = ch_versions.mix(FASTQ_SUBSAMPLE_FQ_SALMON.out.versions)

    FASTQ_SUBSAMPLE_FQ_SALMON
        .out
        .json_info
        .join(ch_strand_fastq.auto_strand)
        .map { meta, json, reads ->
            return [ meta + [ strandedness: getSalmonInferredStrandedness(json) ], reads ]
        }
        .mix(ch_strand_fastq.known_strand)
        .set { ch_strand_inferred_fastq }

    emit:

    reads           = ch_strand_inferred_fastq
    trim_read_count = ch_trim_read_count

    multiqc_files   = ch_multiqc_files.transpose().map{it[1]}
    versions        = ch_versions                     // channel: [ versions.yml ]
}
