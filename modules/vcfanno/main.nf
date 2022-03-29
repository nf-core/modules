process VCFANNO {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::vcfanno=0.3.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcfanno:0.3.3--h9ee0642_0':
        'quay.io/biocontainers/vcfanno:0.3.3--h9ee0642_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta), path(vcf_uncompressed)
    path toml
    path resource_dir

    output:
    tuple val(meta), path("*_annotated.vcf"), emit: vcf
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_vcf = vcf_uncompressed ?: vcf
    """
    ln -sf $resource_dir/* \$(pwd)

    vcfanno \\
        -p $task.cpus \\
        $args \\
        $toml \\
        $input_vcf \\
        > ${prefix}_annotated.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfanno: \$(echo \$(vcfanno 2>&1 | grep version | cut -f3 -d' ' ))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_annotated.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfanno: \$(echo \$(vcfanno 2>&1 | grep version | cut -f3 -d' ' ))
    END_VERSIONS
    """
}
