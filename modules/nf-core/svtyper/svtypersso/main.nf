process SVTYPER_SVTYPERSSO {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svtyper:0.7.1--py_0':
        'biocontainers/svtyper:0.7.1--py_0' }"

    input:
    tuple val(meta), path(bam), path(bam_index), path(vcf)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.vcf") , emit: gt_vcf
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf    = vcf ? "--input_vcf ${vcf}" : ""
    def fasta  = fasta ? "--ref_fasta ${fasta}" : ""
    if ("$vcf" == "${prefix}.vcf") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    svtyper-sso \\
        --bam $bam \\
        $vcf \\
        $fasta \\
        --output_vcf ${prefix}.vcf \\
        --lib_info ${prefix}.json \\
        --cores $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svtyper: \$(echo \$(svtyper-sso -h 2>&1 | grep "version:" | sed 's/^version: v//'))
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.json
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svtyper: \$(echo \$(svtyper-sso -h 2>&1 | grep "version:" | sed 's/^version: v//'))
    END_VERSIONS
    """
}
