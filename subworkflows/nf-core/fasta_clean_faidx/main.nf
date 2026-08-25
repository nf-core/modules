include { GUNZIP                        } from "../../../modules/nf-core/gunzip/main"
include { SEQKIT_SEQ                    } from '../../../modules/nf-core/seqkit/seq/main'
include { SEQKIT_REPLACE as SEQKIT_DOTS } from '../../../modules/nf-core/seqkit/replace/main'
include { SAMTOOLS_FAIDX                } from "../../../modules/nf-core/samtools/faidx/main"
include { SAMTOOLS_DICT                 } from "../../../modules/nf-core/samtools/dict/main"

workflow FASTA_CLEAN_FAIDX {
    take:
    ch_reference          // channel.of( [meta], reference )
    val_get_chromsizes    // boolean: emit chromsizes
    val_get_dict          // boolean: emit dict

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
    //         ALSO RETURN ONLY THE ID OF A SEQUENCE NOT THE FULL HEADER
    //
    SEQKIT_SEQ (
        unzipped_reference
    )


    //
    // MODULE: REPLACE `.` IN HEADERS WITH `_` (THIS CAN SOMETIMES BE IMPORTANT INFORMATION)
    //
    SEQKIT_DOTS (
        SEQKIT_SEQ.out.fastx,
        "fasta"
    )


    //
    // MODULE: GENERATE INDEX OF REFERENCE FASTA
    //         OPTIONALLY EMIT CHROMOSOME SIZES FILE
    //
    SAMTOOLS_FAIDX (
        SEQKIT_DOTS.out.fastx.map { meta, file -> [meta, file, []] },
        val_get_chromsizes
    )

    // LOGIC: GENERATE A SEQUENCE DESCRIPTION FILE FROM THE FAI
    ch_stats_json = SAMTOOLS_FAIDX.out.fai
        .map { meta, fai ->
            def lengths = fai.readLines()
                .findAll { line -> line.trim() }
                .collect { line -> line.split('\t')[1].toLong() }

            def sequence_map = [
                n_sequences : lengths.size(),
                total_length: lengths.sum() ?: 0L,
            ]
            if (lengths) {
                sequence_map.max_length = lengths.max()
            }

            def outfile = file("${params.outdir}/genome_stats/${meta.id}.genome_stats.json")
            outfile.parent.mkdirs()
            outfile.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(sequence_map))


            tuple(meta, outfile)
        }

    //
    // MODULE: GENERATE A SAMTOOLS DICT FILE BASED ON THE CORRECTED FASTA FILE
    //
    SAMTOOLS_DICT (
        SEQKIT_DOTS.out.fastx.filter { meta, file -> val_get_dict }
    )


    emit:
    reference               = SEQKIT_DOTS.out.fastx
    fai                     = SAMTOOLS_FAIDX.out.fai
    sizes                   = SAMTOOLS_FAIDX.out.sizes
    dict                    = SAMTOOLS_DICT.out.dict
    sequence_description    = ch_stats_json
}
