process DECOUPLER_PKN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/saezlab/publish-packages/decoupler:sha-2f65a0d"

    input:
    tuple val(meta), val(resource), val(organism)

    output:
    tuple val(meta), path("*.csv") , emit: network
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix   = task.ext.prefix ?: "${meta.id}_decoupler_network"
    """
    #!/usr/bin/env python3

    import decoupler as dc
    import pandas as pd

    resource_dict = {'dorothea' : dc.get_dorothea,
                    'collectri' : dc.get_collectri,
                    'progeny'   : dc.get_progeny}

    resource_type = "${resource}".lower()
    org = "${organism}".lower()

    get_net_df = resource_dict.get(resource_type)

    net = get_net_df(organism = org)
    net.to_csv("${prefix}.csv")

    ## VERSIONS FILE
    with open('versions.yml', 'a') as version_file:
        version_file.write('"${task.process}":' + "\\n")
        version_file.write("\tdecoupler-py: " + dc.__version__ + "\\n")
    """
    stub:
    prefix   = task.ext.prefix ?: "${meta.id}_decoupler_network"
    """
    #!/usr/bin/env python3
    open("${prefix}.csv", 'a').close()

    with open('versions.yml', 'a') as version_file:
        version_file.write('"${task.process}":' + "\\n")
        version_file.write("\tdecoupler-py: " + dc.__version__ + "\\n")
    """
}
