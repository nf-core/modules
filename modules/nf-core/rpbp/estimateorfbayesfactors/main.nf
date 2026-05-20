process RPBP_ESTIMATEORFBAYESFACTORS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/14/146c3f15abf184a5ec13531d2a040ba7b9235c1091723aa37c7a119817411367/data' :
        'community.wave.seqera.io/library/rpbp:4.0.1--71297b462026e13b' }"

    input:
    tuple val(meta),  path(profiles)
    tuple val(meta4), path(orfs_genomic_bed)
    tuple val(meta2), path(translated_models,   stageAs: 'translated_models/*')
    tuple val(meta3), path(untranslated_models, stageAs: 'untranslated_models/*')

    output:
    tuple val(meta), path("${prefix}.bayes-factors.bed.gz"), emit: bayes_factors
    tuple val("${task.process}"), val('rpbp'), eval('python -c "import rpbp; print(rpbp.__version__)"'), emit: versions_rpbp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def translated_flag   = translated_models   ? "--translated-models ${translated_models.collect { it.toString() }.join(' ')}"     : '--translated-models $TRANSLATED_MODELS'
    def untranslated_flag = untranslated_models ? "--untranslated-models ${untranslated_models.collect { it.toString() }.join(' ')}" : '--untranslated-models $UNTRANSLATED_MODELS'
    def need_bundled      = !translated_models || !untranslated_models
    def bundled_setup     = need_bundled ? '''
    RPBP_MODELS_BASE=$(python3 -c "import os, inspect, rpbp; print(os.path.join(os.path.dirname(inspect.getfile(rpbp)), 'models'))")
    TRANSLATED_MODELS=$(ls $RPBP_MODELS_BASE/translated/*.stan | xargs)
    UNTRANSLATED_MODELS=$(ls $RPBP_MODELS_BASE/untranslated/*.stan | xargs)
    ''' : ''
    """
    ${bundled_setup}
    estimate-orf-bayes-factors \\
        ${profiles} \\
        ${orfs_genomic_bed} \\
        ${prefix}.bayes-factors.bed.gz \\
        ${translated_flag} \\
        ${untranslated_flag} \\
        --num-cpus ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.bayes-factors.bed.gz
    """
}
