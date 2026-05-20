process RPBP_ESTIMATEMETAGENEBAYESFACTORS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/14/146c3f15abf184a5ec13531d2a040ba7b9235c1091723aa37c7a119817411367/data' :
        'community.wave.seqera.io/library/rpbp:4.0.1--71297b462026e13b' }"

    input:
    tuple val(meta),  path(profile_csv)
    tuple val(meta2), path(periodic_models,    stageAs: 'periodic_models/*')
    tuple val(meta3), path(nonperiodic_models, stageAs: 'nonperiodic_models/*')

    output:
    tuple val(meta), path("${prefix}.metagene-periodicity-bayes-factors.csv.gz"), emit: bayes_factors
    tuple val("${task.process}"), val('rpbp'), eval('python -c "import rpbp; print(rpbp.__version__)"'), emit: versions_rpbp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def periodic_paths    = periodic_models    ? periodic_models.join(' ')    : ''
    def nonperiodic_paths = nonperiodic_models ? nonperiodic_models.join(' ') : ''
    """
    PERIODIC="${periodic_paths}"
    NONPERIODIC="${nonperiodic_paths}"
    if [ -z "\$PERIODIC" ] || [ -z "\$NONPERIODIC" ]; then
        RPBP_MODELS_BASE=\$(python3 -c "import os, inspect, rpbp; print(os.path.join(os.path.dirname(inspect.getfile(rpbp)), 'models'))")
        [ -z "\$PERIODIC" ]    && PERIODIC=\$(ls "\$RPBP_MODELS_BASE"/periodic/*.stan)
        [ -z "\$NONPERIODIC" ] && NONPERIODIC=\$(ls "\$RPBP_MODELS_BASE"/nonperiodic/*.stan)
    fi

    estimate-metagene-profile-bayes-factors \\
        ${profile_csv} \\
        ${prefix}.metagene-periodicity-bayes-factors.csv.gz \\
        --periodic-models \$PERIODIC \\
        --nonperiodic-models \$NONPERIODIC \\
        --num-cpus ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.metagene-periodicity-bayes-factors.csv.gz
    """
}
