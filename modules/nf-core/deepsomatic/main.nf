process DEEPSOMATIC {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/google/deepsomatic:1.10.0" // Keep in sync with versions_deepsomatic.
    input:
    tuple val(meta), path(input_normal), path(index_normal), path(input_tumor), path(index_tumor), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)

    output:
    tuple val(meta), path("${prefix}.vcf.gz")      ,  emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi")  ,  emit: vcf_tbi
    tuple val(meta), path("${prefix}.g.vcf.gz")    ,  emit: gvcf, optional: true
    tuple val(meta), path("${prefix}.g.vcf.gz.tbi"),  emit: gvcf_tbi, optional: true
    tuple val("${task.process}"), val("deepsomatic"), val("1.10.0"), emit: versions_deepsomatic, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPSOMATIC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def regions = intervals ? "--regions=${intervals}" : ""

    """
    run_deepsomatic \\
        --ref=${fasta} \\
        --reads_normal=${input_normal} \\
        --reads_tumor=${input_tumor} \\
        --output_vcf=${prefix}.vcf.gz \\
        ${args} \\
        ${regions} \\
        --intermediate_results_dir=tmp \\
        --num_shards=${task.cpus}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPSOMATIC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def gvcf = (args =~ /(?:^|\s)--output_gvcf(?:=|\s)/) ? "echo \"\" | gzip > ${prefix}.g.vcf.gz; touch ${prefix}.g.vcf.gz.tbi" : ""
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    ${gvcf}
    """
}
