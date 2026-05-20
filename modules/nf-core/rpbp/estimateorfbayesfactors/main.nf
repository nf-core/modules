process RPBP_ESTIMATEORFBAYESFACTORS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3a/3a8aa95ce76934f6269b2d8cbdd3d57c13db029c704152975b2315e35b7a2154/data' :
        'community.wave.seqera.io/library/rpbp_star:247a8ae84a6babfb' }"

    input:
    tuple val(meta), path(profiles)
    path  orfs_genomic_bed

    output:
    tuple val(meta), path("${prefix}.bayes-factors.bed.gz"), emit: bayes_factors
    tuple val("${task.process}"), val('rpbp'), eval('python -c "import rpbp; print(rpbp.__version__)"'), emit: versions_rpbp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    RPBP_MODELS_BASE=\$(python3 -c "import os, inspect, rpbp; print(os.path.join(os.path.dirname(inspect.getfile(rpbp)), 'models'))")
    TRANSLATED_MODELS=\$(ls \$RPBP_MODELS_BASE/translated/*.stan | xargs)
    UNTRANSLATED_MODELS=\$(ls \$RPBP_MODELS_BASE/untranslated/*.stan | xargs)

    estimate-orf-bayes-factors \\
        ${profiles} \\
        ${orfs_genomic_bed} \\
        ${prefix}.bayes-factors.bed.gz \\
        --translated-models \$TRANSLATED_MODELS \\
        --untranslated-models \$UNTRANSLATED_MODELS \\
        --num-cpus ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.bayes-factors.bed.gz
    """
}
