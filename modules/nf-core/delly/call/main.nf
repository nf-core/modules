process DELLY_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::delly=1.1.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/delly:1.1.6--ha41ced6_0' :
        'quay.io/biocontainers/delly:1.1.6--ha41ced6_0' }"

    input:
    tuple val(meta), path(input), path(input_index), path(exclude_bed)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.{bcf,vcf}")  , emit: bcf
    tuple val(meta), path("*.csi")      , emit: csi, optional:true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "bcf"
    def bcf_output = suffix == "bcf" ? "--outfile ${prefix}.bcf" : ""
    def vcf_output = suffix == "vcf" ? "> ${prefix}.vcf" : ""
    def exclude = exclude_bed ? "--exclude ${exclude_bed}" : ""
    """
    delly \\
        call \\
        ${args} \\
        ${bcf_output} \\
        --genome ${fasta} \\
        ${exclude} \\
        ${input} \\
        ${vcf_output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$( echo \$(delly --version 2>&1) | sed 's/^.*Delly version: v//; s/ using.*\$//')
    END_VERSIONS
    """
}
