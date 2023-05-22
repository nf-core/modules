process SVTYPER_SVTYPER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::svtyper=0.7.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svtyper:0.7.0--py27h24bf2e0_1':
        'biocontainers/svtyper:0.7.0--py27h24bf2e0_1' }"

    input:
    tuple val(meta), path(bam), path(bam_index), path(vcf)
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.vcf") , emit: gt_vcf
    tuple val(meta), path("*.bam") , emit: bam
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf  = vcf ? "--input_vcf ${vcf}" : ""
    if ("$vcf" == "${prefix}.vcf") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    svtyper \\
        $vcf \\
        --bam $bam \\
        --lib_info ${prefix}.json \\
        --output_vcf ${prefix}.vcf \\
        --ref_fasta $fasta \\
        --write_alignment ${prefix}.bam \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svtyper: \$(echo \$(svtyper -h 2>&1 | grep "version:" | sed 's/^version: v//'))
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.json
    touch ${prefix}.vcf
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svtyper: \$(echo \$(svtyper -h 2>&1 | grep "version:" | sed 's/^version: v//'))
    END_VERSIONS
    """
}
