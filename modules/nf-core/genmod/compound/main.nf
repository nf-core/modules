process GENMOD_COMPOUND {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::genmod=3.7.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genmod:3.7.4--pyh5e36f6f_0':
        'biocontainers/genmod:3.7.4--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(input_vcf)

    output:
    tuple val(meta), path("*_compound.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    genmod \\
        compound \\
        $args \\
        --processes ${task.cpus} \\
        --outfile ${prefix}_compound.vcf \\
        $input_vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmod: \$(echo \$(genmod --version 2>&1) | sed 's/^.*genmod version: //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_compound.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmod: \$(echo \$(genmod --version 2>&1) | sed 's/^.*genmod version: //' )
    END_VERSIONS
    """
}
