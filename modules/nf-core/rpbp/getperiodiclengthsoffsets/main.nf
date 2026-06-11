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
    tuple val(meta), path("${prefix}.lengths-offsets.tsv"), emit: lengths_offsets
    path "versions.yml"                                            , emit: versions_rpbp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Periodic-length filter thresholds, space-separated:
    //   <min_metagene_profile_count> <min_metagene_bf_mean> <max_metagene_bf_var> <min_metagene_bf_likelihood>
    // Defaults mirror rpbp.defaults.metagene_options. Any token may be "None"
    // to defer to upstream's per-slot default-filter handling.
    def filter_args = (task.ext.args ?: '1000 5 None 0.5').tokenize(' ')
    prefix = task.ext.prefix ?: "${meta.id}"
    min_count   = filter_args[0]
    min_bf_mean = filter_args[1]
    max_bf_var  = filter_args[2]
    min_bf_lik  = filter_args[3]
    template 'get_periodic_lengths_and_offsets.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.lengths-offsets.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
        rpbp: \$(python -c "import rpbp; print(rpbp.__version__)")
    END_VERSIONS
    """
}
