// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { FASTQ_FETCH_CLEAN_HOSTILE } from '../fastq_fetch_clean_hostile/main'
include { FASTQ_INDEX_FILTER_DEACON } from '../fastq_index_filter_deacon/main'

workflow FASTQ_DECONTAMINATE_DEACON_HOSTILE {

    take:
    // TODO nf-core: edit input (take) channels
    ch_reads                // channel: [ val(meta), [ reads ] ]
    ch_fasta                // channel: [ val(meta), [ fasta ] ]
    ch_reference            // channel: [ val(reference_name), path(reference_dir) ] (optional)
    index_name              // val (optional)
    decontaminator          // string (enum): 'hostile' or 'deacon'

    main:

    ch_versions = Channel.empty()

    if (decontaminator != "hostile" && decontaminator != "deacon"){
        error("Unknown decontaminator '${decontaminator}'")
    }

    // Fastq decontamination
    if (decontaminator == "hostile") {
        FASTQ_FETCH_CLEAN_HOSTILE (
            ch_reads,
            ch_reference,
            index_name
        )
        fastq_filtered = FASTQ_FETCH_CLEAN_HOSTILE.out.fastq
        ch_versions = ch_versions.mix(FASTQ_FETCH_CLEAN_HOSTILE.out.versions)
        additional_files = FASTQ_FETCH_CLEAN_HOSTILE.out.reference
            .concat(FASTQ_FETCH_CLEAN_HOSTILE.out.json)

    } else if (decontaminator == "deacon") {
        FASTQ_INDEX_FILTER_DEACON (
            ch_fasta.join(ch_reads)
        )
        fastq_filtered = FASTQ_INDEX_FILTER_DEACON.out.fastq_filtered
        ch_versions = ch_versions.mix(FASTQ_INDEX_FILTER_DEACON.out.versions)
        additional_files = FASTQ_INDEX_FILTER_DEACON.out.index
            .join(FASTQ_INDEX_FILTER_DEACON.out.summary)
    }


    emit:
    fastq_filtered      = fastq_filtered    // channel: [ val(meta), [ fastq ] ]
    additional_files    = additional_files  // channel: file
    versions            = ch_versions       // channel: [ versions.yml ]
}
