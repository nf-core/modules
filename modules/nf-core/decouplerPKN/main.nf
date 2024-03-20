process GET_DECOUPLER_PKN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-decoupler:2.8.0--r43hdfd78af_0':
        'biocontainers/bioconductor-decoupler:2.8.0--r43hdfd78af_0'}"

    input:
    
    tuple val(meta) ,val(resource), val(organism)

    output:
    
    tuple val(meta), path("decoupler_network.csv")          , emit: network
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    
    """
    #!/usr/bin/env python3
    
    import decoupler as dc
    import pandas as pd

    resource_dict = {'dorothea' : dc.get_dorothea, 'collectri' : dc.get_collectri,
    'progeny' : dc.get_progeny}

    resource_type = "${resource}".lower()
    org = "${organism}".lower()

    get_net_df = resource_dict.get(resource_type)
    
    net = get_net_df(organism = org)
    net.to_csv('decoupler_network.csv')

    ## VERSIONS FILE
    with open('versions.yml', 'a') as version_file:
        version_file.write('"${task.process}":' + "\\n")
        version_file.write("\tdecoupler-py: " + dc.__version__ + "\\n")
    """
}
