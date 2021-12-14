process VCFANNO_LUA {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::vcfanno=0.3.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcfanno:0.3.3--h9ee0642_0':
        'quay.io/biocontainers/vcfanno:0.3.3--h9ee0642_0' }"

    input:
    tuple val(meta), path(vcf)
    path vcfanno_config
    path vcfanno_functions

    output:
    tuple val(meta), path("*.vcf"), emit: vcf

    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def lua_functions = vcfanno_functions ?: ""
    """
    vcfanno \\
        -p $task.cpus \\
        $args \\
        -lua \\
        ${vcfanno_config} \\
        $lua_functions \\
        ${vcf} \\
        > ${prefix}_annotated.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfanno: \$(echo \$(vcfanno 2>&1) | grep version | cut -f3 -d' ' ))
    END_VERSIONS
    """
}
