process PARABRICKS_DBSNP {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'
    stageInMode 'copy' // needed by the module to work properly - might be removed when this is fixed upstream

    container "nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1"

    input:
    tuple val(meta), path(vcf_file), path(dbsnp_file), path(tabix_file)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def num_gpus = task.accelerator ? "--num-gpus $task.accelerator.request" : ''
    """
    pbrun \\
        dbsnp \\
        --in-vcf $vcf_file \\
        --in-dbsnp-file $dbsnp_file \\
        --out-vcf ${prefix}.vcf \\
        $num_gpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
