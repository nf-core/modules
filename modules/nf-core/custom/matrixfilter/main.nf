process CUSTOM_MATRIXFILTER {
    tag "$meta"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/48/483e9d9b3b07e5658792d579e230ad40ed18daf7b9ebfb4323c08570f92fd1d5/data'
        : 'community.wave.seqera.io/library/r-base:4.2.1--b0b5476e2e7a0872'}"

    input:
    tuple val(meta), path(abundance)
    tuple val(samplesheet_meta), path(samplesheet)

    output:
    tuple val(meta), path("*.filtered.tsv")             , emit: filtered
    tuple val(meta), path("*.tests.tsv")                , emit: tests
    tuple val(meta), path("*R_sessionInfo.log")         , emit: session_info
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Note: params are passed at line 100 of the template like:
    //
    // opt <- parse_args('$task.ext.args', opt)
    //
    // (new variables defined here don't seem to be available in templates, so
    // we have to access $task directly)
    template 'matrixfilter.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered.tsv
    touch ${prefix}.tests.tsv
    touch ${prefix}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
