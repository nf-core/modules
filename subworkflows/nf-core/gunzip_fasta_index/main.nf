include { GUNZIP                        } from "../../../modules/nf-core/gunzip/main"
include { GAWK as GAWK_UPPER_SEQUENCE   } from '../../../modules/nf-core/gawk/main'
include { SAMTOOLS_FAIDX                } from "../../../modules/nf-core/samtools/faidx/main"

workflow SETUP_FASTA {
    take:
    ch_reference          // channel.of( [meta], reference )
    val_get_chromsizes    // boolean: emit chromsizes

    main:

    //
    // LOGIC: SPLIT THE INPUT REFERENCES INTO THOSE THAT ARE ZIPPED OR UNZIPPED
    //        THIS ENSURES THAT ALL INPUT ARE UNZIPPED FOR DOWNSTREAM PROCESSING
    //
    ch_input = ch_reference
            .branch { _meta, file ->
                zipped: file.name.endsWith('.gz')
                unzipped: !file.name.endsWith('.gz')
            }


    //
    // MODULE: UNZIP INPUTS IF NEEDED
    //
    GUNZIP (
        ch_input.zipped
    )


    //
    // LOGIC: MIX CHANELS WHICH MAY OR MAY NOT BE EMPTY INTO A SINGLE QUEUE CHANNEL
    //
    unzipped_reference = ch_input.unzipped
        .mix(GUNZIP.out.gunzip)


    //
    // MODULE: UPPERCASE THE REFERENCE SEQUENCE
    //         AWK = IF HEADER PASS ELSE CONVERT LINE TO UPPER CASE
    //
    ch_upper_sequence = channel.of('''\
        /^>/ { print; next }
        !/^>/ { print toupper(\$0) }
        '''.stripIndent()
    ).collectFile(name: "uppercase_sequence.awk", cache: true)
    .collect()

    GAWK_UPPER_SEQUENCE(
        unzipped_reference,
        ch_upper_sequence,
        false,
    )


    //
    // MODULE: GENERATE INDEX OF REFERENCE FASTA
    //         OPTIONALLY EMIT CHROMOSOME SIZES FILE
    //
    SAMTOOLS_FAIDX (
        GAWK_UPPER_SEQUENCE.out.output.map { meta, file -> [meta, file, []] },
        val_get_chromsizes
    )


    emit:
    reference = GAWK_UPPER_SEQUENCE.out.output
    fai       = SAMTOOLS_FAIDX.out.fai
    sizes     = SAMTOOLS_FAIDX.out.sizes
}
