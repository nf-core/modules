process DECOUPLER {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cf/cf3578c5f4a5ad40d10b535d496daa29f47e908fbbbc789f85e3341b8f9341ba/data'
        : 'community.wave.seqera.io/library/decoupler-py:1.8.0--b139c6574c0ac2a6'}"

    input:
    tuple val(meta), path(mat)
    path net

    output:
    tuple val(meta), path("*estimate__decoupler.tsv"), emit: dc_estimate
    tuple val(meta), path("*pvals__decoupler.tsv"), emit: dc_pvals
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: "{}"
    """
    #!/usr/bin/env python3
    import os
    import pandas as pd

    os.environ["NUMBA_CACHE_DIR"] = "./tmp"

    import decoupler as dc

    methods = ['aucell', 'gsea', 'gsva', 'mdt', 'mlm', 'ora', 'udt',
        'ulm', 'viper', 'wmean', 'wsum']

    mat = pd.read_csv("${mat}", sep="\t", index_col=0)
    net = pd.read_csv("${net}", sep="\t", index_col=0)

    # Parsing arguments
    args = ${args}
    parsedargs = {'args': {}}

    for k, v in args.items():
        # Specific method argument
        if k.split('_')[0] in methods:
            meth = k.split('_')[0]
            arg = '_'.join(k.split('_')[1:])

            if meth not in args['args'].keys():
                parsedargs['args'][meth] = {arg: v}
            else:
                parsedargs['args'][meth].update({arg: v})

        # Generic argument
        else:
            parsedargs[k] = v


    results = dc.decouple(
        mat=mat,
        net=net,
        **parsedargs
    )

    for result in results:
        results[result].to_csv(result + "__decoupler.tsv", sep="\t")

    ## VERSIONS FILE
    with open('versions.yml', 'a') as version_file:
        version_file.write('"${task.process}":' + "\\n")
        version_file.write("\tdecoupler-py: " + dc.__version__ + "\\n")
    """
}
