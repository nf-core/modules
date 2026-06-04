process HOSTILE_FETCH {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7c/7caca3a47606de8e3460b35823193a471272aa6ab7cfafbf9aabf4615c9fa181/data'
        : 'community.wave.seqera.io/library/hostile:2.0.2--a7f5e5d341b6b94b'}"

    input:
    val index_name

    output:
    tuple val(index_name), path('reference/'), emit: reference
    tuple val("${task.process}"), val('hostile'), eval("hostile --version"), emit: versions_hostile, topic: versions


    script:
    def args = task.ext.args ?: ''
    """
    mkdir reference/
    export HOSTILE_CACHE_DIR=./reference

    hostile \\
        index \\
        fetch \\
        --name ${index_name} \\
        ${args}
    """

    stub:
    """
    mkdir reference/
    export HOSTILE_CACHE_DIR=./reference
    touch reference/human-t2t-hla.1.bt2
    touch reference/human-t2t-hla.2.bt2
    touch reference/human-t2t-hla.3.bt2
    touch reference/human-t2t-hla.4.bt2
    touch reference/human-t2t-hla.rev.1.bt2
    touch reference/human-t2t-hla.rev.2.bt2
    touch reference/human-t2t-hla.mmi
    """
}
