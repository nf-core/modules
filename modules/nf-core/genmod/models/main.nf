process GENMOD_MODELS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genmod:3.9--pyhdfd78af_0':
        'biocontainers/genmod:3.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input_vcf), path (fam)
    path (reduced_penetrance)

    output:
    tuple val(meta), path("*_models.vcf"), emit: vcf
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def family_file =  fam ? "--family_file ${fam}" : ""
    def pen_file    = reduced_penetrance ? "--reduced_penetrance ${reduced_penetrance}" : ""
    """
    genmod \\
        models \\
        $args \\
        $pen_file \\
        $family_file \\
        --processes ${task.cpus} \\
        --outfile ${prefix}_models.vcf \\
        $input_vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmod: \$(echo \$(genmod --version 2>&1) | sed 's/^.*genmod version: //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_models.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmod: \$(echo \$(genmod --version 2>&1) | sed 's/^.*genmod version: //' )
    END_VERSIONS
    """
}
