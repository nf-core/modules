process PARABRICKS_DBSNP {
    tag "$meta.id"
    label 'process_high'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.3.0-1"

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
    def ann_vcf = "${prefix}.vcf"
    """

    pbrun \\
        dbsnp \\
        --in-vcf $vcf_file \\
        --in-dbsnp-file $dbsnp_file \\
        --out-vcf $ann_vcf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ann_vcf = "${prefix}.vcf"
    """
    touch $ann_vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
