process RPBP_ESTIMATEMETAGENEBAYESFACTORS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3a/3a8aa95ce76934f6269b2d8cbdd3d57c13db029c704152975b2315e35b7a2154/data' :
        'community.wave.seqera.io/library/rpbp_star:247a8ae84a6babfb' }"

    input:
    tuple val(meta), path(profile_csv)

    output:
    tuple val(meta), path("${prefix}.metagene-periodicity-bayes-factors.csv.gz"), emit: bayes_factors
    tuple val("${task.process}"), val('rpbp'), eval('python -c "import rpbp; print(rpbp.__version__)"'), emit: versions_rpbp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    RPBP_MODELS_BASE=\$(python3 -c "import os, inspect, rpbp; print(os.path.join(os.path.dirname(inspect.getfile(rpbp)), 'models'))")
    PERIODIC_MODELS=\$(ls \$RPBP_MODELS_BASE/periodic/*.stan | xargs)
    NONPERIODIC_MODELS=\$(ls \$RPBP_MODELS_BASE/nonperiodic/*.stan | xargs)

    estimate-metagene-profile-bayes-factors \\
        ${profile_csv} \\
        ${prefix}.metagene-periodicity-bayes-factors.csv.gz \\
        --periodic-models \$PERIODIC_MODELS \\
        --nonperiodic-models \$NONPERIODIC_MODELS \\
        --num-cpus ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.metagene-periodicity-bayes-factors.csv.gz
    """
}
