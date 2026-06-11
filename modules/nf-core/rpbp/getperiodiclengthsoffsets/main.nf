process RPBP_GETPERIODICLENGTHSOFFSETS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/14/146c3f15abf184a5ec13531d2a040ba7b9235c1091723aa37c7a119817411367/data' :
        'community.wave.seqera.io/library/rpbp:4.0.1--71297b462026e13b' }"

    input:
    tuple val(meta), path(periodic_offsets)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: lengths_offsets
    path "versions.yml"                                            , emit: versions_rpbp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    task_ext_args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.lengths-offsets"
    template 'get_periodic_lengths_and_offsets.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.lengths-offsets"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
        rpbp: \$(python -c "import rpbp; print(rpbp.__version__)")
    END_VERSIONS
    """
}
