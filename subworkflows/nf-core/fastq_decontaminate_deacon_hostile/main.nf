include { FASTQ_FETCH_CLEAN_HOSTILE } from '../fastq_fetch_clean_hostile/main'
include { FASTQ_INDEX_FILTER_DEACON } from '../fastq_index_filter_deacon/main'

workflow FASTQ_DECONTAMINATE_DEACON_HOSTILE {

    take:
    ch_reads                // channel: [ val(meta), [ reads ] ]
    ch_fasta                // channel: [ val(meta), [ fasta ] ] (optional)
    ch_reference            // channel: [ val(reference_name), path(reference_dir) ] (optional)
    index_name              // val (optional)
    decontaminator          // string (enum): 'hostile' or 'deacon'

    main:

    reference = channel.empty()
    json = channel.empty()
    index = channel.empty()
    summary = channel.empty()

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
        reference = FASTQ_FETCH_CLEAN_HOSTILE.out.reference
        json = FASTQ_FETCH_CLEAN_HOSTILE.out.json
    } else if (decontaminator == "deacon") {
        FASTQ_INDEX_FILTER_DEACON (
            ch_fasta.join(ch_reads)
        )
        fastq_filtered = FASTQ_INDEX_FILTER_DEACON.out.fastq_filtered
        index = FASTQ_INDEX_FILTER_DEACON.out.index
        summary = FASTQ_INDEX_FILTER_DEACON.out.summary
    }


    emit:
    fastq_filtered      = fastq_filtered    // channel: [ val(meta), [ fastq ] ]
    reference           = reference         // channel: [ val(reference_name), path(reference_dir) ] (hostile only)
    json                = json              // channel: [ val(meta), [ *.json ] ] (hostile only)
    index               = index             // channel: [ val(meta), [ index ] ] (deacon only)
    summary             = summary           // channel: [ val(meta), [ log ] ] (deacon only)
}
