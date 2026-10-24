process BUSCO_GENERATEPLOT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dc1935be05e41ef9e53df8820a4862bc2f35637af83fbe1dc402fc109824e5a2/data'
        : 'community.wave.seqera.io/library/busco_sepp:28bfd2435aed3cc4' }"

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
