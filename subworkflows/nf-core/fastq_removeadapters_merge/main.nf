// ok for single end adapter removal
include { TRIMMOMATIC   } from '../../../modules/nf-core/trimmomatic/main'
include { CUTADAPT      } from '../../../modules/nf-core/cutadapt/main'
include { TRIMGALORE    } from '../../../modules/nf-core/trimgalore/main'
include { BBMAP_BBDUK   } from '../../../modules/nf-core/bbmap/bbduk/main'

// allows merging of paired end reads, but will work for single end reads as well
include { FASTP         } from '../../../modules/nf-core/fastp/main'
include { ADAPTERREMOVAL} from '../../../modules/nf-core/adapterremoval/main'
include { LEEHOM        } from '../../../modules/nf-core/leehom/main'

// requirs paired end becouse of merging
include { NGMERGE       } from '../../../modules/nf-core/ngmerge/main'

workflow FASTQ_REMOVEADAPTERS_MERGE {

    take:
    reads                    // (meta, reads) - SE or PE
    skip_trimmomatic         // boolean
    skip_cutadapt
    skip_trimgalore
    skip_bbduk
    skip_fastp
    skip_adapterremoval
    skip_leehom
    skip_ngmerge
    do_merge
    adapters_contaminants
    main:

    ch_current = reads
    ch_versions = Channel.empty()

    if (!skip_trimmomatic) {
   //     TRIMMOMATIC.ext.adapters = adapters_contaminants
        TRIMMOMATIC( ch_current )
        ch_current = TRIMMOMATIC.out.trimmed_reads.map { meta, r -> [meta, r] }
        ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)
    }

    if (!skip_cutadapt) {
        CUTADAPT( ch_current )
        ch_current = CUTADAPT.out.reads.map { meta, r -> [meta, r] }
        ch_versions = ch_versions.mix(CUTADAPT.out.versions)
    }

    if (!skip_trimgalore) {
        TRIMGALORE( ch_current )
        ch_current = TRIMGALORE.out.reads.map { meta, r -> [meta, r] }
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
    }

    if (!skip_bbduk) {
        ch_current.view { "DEBUG: BEFORE BBMAP_BBDUK → $it" }
        BBMAP_BBDUK( ch_current, adapters_contaminants )
        ch_current = BBMAP_BBDUK.out.reads.map { meta, r -> [meta, r] }
        ch_versions = ch_versions.mix(BBMAP_BBDUK.out.versions)
        ch_current.view { "DEBUG: AFTER BBMAP_BBDUK → $it" }
    }

    if (!skip_fastp && !do_merge) {
        ch_current.view { "DEBUG: BEFORE FASTP → $it" }
        FASTP( ch_current.map { meta, r ->  tuple(meta, r, adapters_contaminants)},
            false,
            true,
            do_merge)
        ch_current = FASTP.out.reads.map { meta, r -> [meta, r] }
        ch_versions = ch_versions.mix(FASTP.out.versions)
        ch_current.view { "DEBUG: AFTER FASTP → $it" }
    }

    if (!skip_adapterremoval && !do_merge) {
    ADAPTERREMOVAL( ch_current, adapters_contaminants )
    ch_current.view { "DEBUG: BEFORE ADAPTERREMOVAL → $it" }
         ch_current = ADAPTERREMOVAL.out.paired_truncated
    .mix(ADAPTERREMOVAL.out.singles_truncated)
    .map { tup ->
        def meta  = tup[0]
        def r = tup[1]
        [meta, r]
    }
         ch_versions = ch_versions.mix(ADAPTERREMOVAL.out.versions)
         ch_current.view { "DEBUG: AFTER ADAPTERREMOVAL → $it" }
    }

    if (!skip_leehom && !do_merge) {
        ch_current.view { "DEBUG: BEFORE LEEHOM → $it" }
        ch_leehom_input = ch_current.map { meta, r ->
        if (meta.single_end)
            return [meta, r]
        else
            return [meta, [reads[0], reads[1]]]
    }

    LEEHOM(ch_leehom_input)

        ch_leehom_pe =
    LEEHOM.out.unmerged_r1_fq_pass
        .combine(LEEHOM.out.unmerged_r2_fq_pass)
        .map { left, right ->
            def meta = left[0]
            def r1   = left[1]
            def r2   = right[1]
            [meta, [r1, r2]]
        }
        .filter { meta, r -> !meta.single_end }  // only PE

// ---- SE OUTPUT ----
ch_leehom_se =
    LEEHOM.out.fq_pass
        .filter { meta, r -> meta.single_end }
        .map { meta, r -> [meta, r] }

// ---- MERGE CHANNELS ----
ch_current = ch_leehom_pe.mix(ch_leehom_se)

        ch_versions = ch_versions.mix(LEEHOM.out.versions)
    }

    emit:

    trimmed_reads = ch_current               // channel: [ meta, trimmed_reads ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
