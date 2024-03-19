#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DECOUPLER } from '../../../../modules/nf-core/decoupler/main.nf'

process GET_DATA {

    script:
    """
    #!/usr/bin/env python3

    import decoupler as dc

    mat, net = dc.get_toy_data()

    mat.to_csv('mat.csv')
    net.to_csv('net.csv')
    """
}

workflow {
    GET_DATA()

    //DECOUPLER ( input )
}
