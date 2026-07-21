process BUSCO_GENERATEPLOT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c6/c607f319867d96a38c8502f751458aa78bbd18fe4c7c4fa6b9d8350e6ba11ebe/data'
        : 'community.wave.seqera.io/library/busco_sepp:f2dbc18a2f7a5b64'}"

    input:
    path short_summary_txt, stageAs: 'busco/*'

    output:
    path '*.png'        , emit: png
    tuple val("${task.process}"), val('busco'), eval("busco --version 2> /dev/null | sed 's/BUSCO //g'"), emit: versions_busco, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/busco/plot instead.

Reason:
The standalone generate_plot.py script was removed in BUSCO v6.0.0. Plotting is now
built into BUSCO via 'busco --plot', implemented by the busco/plot module.
"""
    assert false: deprecation_message
}
