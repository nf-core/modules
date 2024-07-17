process DEEPVARIANT_MAKEEXAMPLES {
    tag "$meta.id"
    label 'process_high'

    //Conda is not supported at the moment
    container "nf-core/deepvariant:1.5.0"

    input:
    tuple val(meta), path(input), path(index), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)

    output:
    // Output a unique ID (hash) for this task, that can be used to join the gvcf and variant calls later.
    tuple val("${task.hash}"), val(meta), val("${prefix}.examples.tfrecord@${task.cpus}.gz"), path("${prefix}.examples.tfrecord-*-of-*.gz"),    emit: examples
    tuple val("${task.hash}"), val(meta), val("${prefix}.gvcf.tfrecord@${task.cpus}.gz"), path("${prefix}.gvcf.tfrecord-*-of-*.gz"),            emit: gvcf
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def regions = intervals ? "--regions=${intervals}" : ""

    """
    seq 0 ${task.cpus - 1} | parallel -q --halt 2 --line-buffer /opt/deepvariant/bin/make_examples \\
        --mode=calling \\
        --ref=${fasta} \\
        --reads=${input} \\
        --examples "./${prefix}.examples.tfrecord@${task.cpus}.gz" \\
        --channels "insert_size" \\
        --gvcf "./${prefix}.gvcf.tfrecord@${task.cpus}.gz" \\
        ${regions} \\
        ${args} \\
        --task {}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant_makeexamples: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """
}
