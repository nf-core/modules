process METASPACE_SUBMIT {
    label 'process_low'
    container 'docker.io/bwadie/metaspace_converter:latest'

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"

    input:
    path imzml
    path ibd
    path config

    output:
    path("ds_id.txt"), emit: ds_id
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    template 'submit.py'

    stub:

    """
    touch ds_id.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
        metaspace: 2.0.9
        yaml: 6.0.2
    END_VERSIONS
    """
}
