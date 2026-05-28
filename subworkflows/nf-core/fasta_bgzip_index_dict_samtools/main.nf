include { HTSLIB_BGZIPTABIX  } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { SAMTOOLS_DICT      } from '../../../modules/nf-core/samtools/dict/main'
include { SAMTOOLS_FAIDX     } from '../../../modules/nf-core/samtools/faidx/main'

workflow FASTA_BGZIP_INDEX_DICT_SAMTOOLS {

    take:
    ch_fasta // channel: [ val(meta), fasta ]

    main:

    HTSLIB_BGZIPTABIX (
        ch_fasta.map { meta, fasta -> [meta, fasta, [], []] },
        'compress',
        [],
        []
    )

    SAMTOOLS_FAIDX (
        HTSLIB_BGZIPTABIX.out.output.map {meta, fasta -> [meta, fasta, []]},
        true
    )

    SAMTOOLS_DICT (
        HTSLIB_BGZIPTABIX.out.output
    )

    ch_joined = HTSLIB_BGZIPTABIX.out.output
        .join(SAMTOOLS_FAIDX.out.fai)
        .join(SAMTOOLS_FAIDX.out.gzi)
        .join(SAMTOOLS_FAIDX.out.sizes)
        .join(SAMTOOLS_DICT.out.dict)

    emit:
    fasta_fai_gzi_dict = ch_joined             // channel: [ val(meta),  fasta.gz, fai, gzi, sizes, dict ]
}
