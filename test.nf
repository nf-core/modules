include { PHOLD_RUN } from './modules/nf-core/phold/run/main'

workflow {

    ch_input = Channel.of(
        [
            [ id: 'NC_043029' ],
            file('../phold/tests/test_data/NC_043029.gbk')
        ]
    )

    PHOLD_RUN(ch_input)
}