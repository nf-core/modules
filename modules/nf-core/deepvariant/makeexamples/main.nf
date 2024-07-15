process DEEPVARIANT_MAKE_EXAMPLES {
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
    tuple val(meta), val(intervals), val("examples.tfrecord@${task.cpus}.gz"), path("examples.tfrecord-*-of-*.gz"), emit: examples
    tuple val(meta), val(intervals), val("gvcf.tfrecord@${task.cpus}.gz"), path("gvcf.tfrecord-*-of-*.gz"),      emit: gvcf
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def regions = intervals ? "--regions=${intervals}" : ""

    """
    seq 0 ${task.cpus - 1} | parallel -q --halt 2 --line-buffer /opt/deepvariant/bin/make_examples \\
        --mode=calling \\
        --ref=${fasta} \\
        --reads=${input} \\
        --examples "./examples.tfrecord@${task.cpus}.gz" \\
        --channels "insert_size" \\
        --gvcf "./gvcf.tfrecord@${task.cpus}.gz" \\
        ${regions} \\
        ${args} \\
        --task {}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant_makeexamples: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """
}
