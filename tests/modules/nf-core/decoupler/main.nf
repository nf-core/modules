#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DECOUPLER } from '../../../../modules/nf-core/decoupler/main.nf'

process GET_DATA {
    conda "conda-forge::decoupler-py=1.6.0"
    // container = "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'ghcr.io/saezlab/publish-packages/decoupler:sha-2f65a0d' : ''}"

    output:
    path('mat.csv'), emit: mat
    path('net.csv'), emit: net

    script:
    """
    #!/usr/bin/env python3

    import decoupler as dc

    mat, net = dc.get_toy_data()

    mat.to_csv('mat.csv', sep='\t')
    net.to_csv('net.csv', sep='\t')
    """
}

workflow {
    GET_DATA()


    Channel
        .value('ulm')
        .set{method}

    DECOUPLER (Channel.value('toy'), GET_DATA.out.mat, GET_DATA.out.net, method)
}
