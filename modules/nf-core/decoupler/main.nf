process DECOUPLER {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::decoupler-py=1.6.0"
    container = "ghcr.io/saezlab/publish-packages/decoupler:sha-2f65a0d"

    input:
    tuple val(meta), path(mat)
    path(net)

    output:
    tuple val(meta), path("*estimate__decoupler.tsv"), emit: dc_estimate
    tuple val(meta), path("*pvals__decoupler.tsv"), emit: dc_pvals
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '{}'
    def methods = task.ext.methods ?: 'None'
    def source = task.ext.source ?: "source"
    def target = task.ext.target ?: "target"
    def weight = task.ext.weight ?: "weight"
    def min_n = task.ext.min_n ?: 5
    def dense = task.ext.verbose ?: 'False'
    def consensus = task.ext.verbose ?: 'False'
    def verbose = task.ext.verbose ?: 'False'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    import decoupler as dc
    import pandas as pd

    mat = pd.read_csv("${mat}", sep="\t", index_col=0)
    net = pd.read_csv("${net}", sep="\t", index_col=0)

    results = dc.decouple(
        mat=mat,
        net=net,
        methods=${methods},
        source="${source}",
        target="${target}",
        weight="${weight}",
        min_n=int(${min_n}),
        dense=${dense},
        consensus=${consensus},
        verbose=${verbose},
        args=${args}
    )

    for result in results:
        results[result].to_csv(result + "__decoupler.tsv", sep="\t")

    ## VERSIONS FILE
    with open('versions.yml', 'a') as version_file:
        version_file.write('"${task.process}":' + "\\n")
        version_file.write("\tdecoupler-py: " + dc.__version__ + "\\n")
    """
}
