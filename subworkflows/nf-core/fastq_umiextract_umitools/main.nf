include { CAT_FASTQ }        from '../../../modules/nf-core/cat/fastq/main'
include { SEQKIT_SPLIT2 }    from '../../../modules/nf-core/seqkit/split2/main'
include { UMITOOLS_EXTRACT } from '../../../modules/nf-core/umitools/extract/main'

workflow FASTQ_UMIEXTRACT_UMITOOLS {

    take:
    ch_fastqs // channel: [ val(meta), [ fastqs ] ]

    main:

    ch_versions = Channel.empty()

    SEQKIT_SPLIT2(ch_fastqs)

    SEQKIT_SPLIT2
    .out
    .reads
    .map {meta, reads -> [meta, reads.flatten()]}
    .transpose()
    .map{ meta, reads ->
        clone = meta.clone()
        clone['part'] = reads.name.replaceAll("_R?([12]).part","part")
        // this is joining the R1s and R2s from the same part together.
        // Does this need to be smarter, currently requiring _1/_2 or _R1/_R2 in the name

        [ clone, reads ]
    }
    .groupTuple( by : [0])
    .set { ch_fastqs }

    UMITOOLS_EXTRACT ( ch_fastqs )
    ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions)

    UMITOOLS_EXTRACT
    .out
    .reads
    .map{ meta, reads ->
        clone = meta.clone()
        clone['part'] = ""
        [ clone, reads ]
    }
    .groupTuple(by: [0])
    .map { meta, reads -> [ meta, reads.flatten() ] }
    .set { ch_fastqs }

    CAT_FASTQ(ch_fastqs)
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    emit:

    fastqs   = CAT_FASTQ.out.reads             // channel: [ val(meta), [ fastqs ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}

