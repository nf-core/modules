process GRAPHTYPER_VCFCONCATENATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/graphtyper:2.7.7--h7594796_1':
        'biocontainers/graphtyper:2.7.7--h7594796_1' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}.vcf.gz")    , emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("$vcf" == "${prefix}.vcf.gz") {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }
    input_vcfs = vcf.collate(1000).collect{it.join(' ')} // Batching needed because if there are too many VCFs the shell cannot run the command
    commands = input_vcfs.withIndex().collect{ batch_vcfs, index ->
        "graphtyper vcf_concatenate ${batch_vcfs} ${args} --output=${prefix}_subset_${index}.vcf.gz"
    }.join('\n    ')
    """
    # Run each batch of VCFs
    ${commands}

    # Combine the VCF for each batch into a single output VCF
    graphtyper vcf_concatenate \\
        ${prefix}_subset_*.vcf.gz \\
        $args \\
        --write_tbi \\
        --output=${prefix}.vcf.gz

    # Delete batch VCFs
    rm ${prefix}_subset_*.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graphtyper: \$(graphtyper --help | tail -n 1 | sed 's/^   //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("$vcf" == "${prefix}.vcf.gz") {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }
    """
    echo | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graphtyper: \$(graphtyper --help | tail -n 1 | sed 's/^   //')
    END_VERSIONS
    """
}
