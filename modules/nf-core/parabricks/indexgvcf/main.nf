process PARABRICKS_INDEXGVCF {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1"

    input:
    tuple val(meta), path(gvcf)

    output:
    tuple val(meta), path("*.{idx,tbi}") , emit: gvcf_index
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def num_gpus = task.accelerator ? "--num-gpus $task.accelerator.request" : ''
    """
    pbrun \\
        indexgvcf \\
        --input $gvcf \\
        $num_gpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_cmd = gvcf.any{ it.name.endsWith(".gz") } ? "touch ${prefix}.g.vcf.gz.tbi" : "touch ${prefix}.g.vcf.idx"
    """
    $output_cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
