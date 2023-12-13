process TCOFFEE_IRMSD {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/t-coffee:13.46.0.919e8c6b--hfc96bf3_0':
        'biocontainers/t-coffee:13.46.0.919e8c6b--hfc96bf3_0' }"

    input:
    tuple  val(meta),  file (msa)
    tuple  val(meta2), file(template), file(structures)

    output:
    tuple val(meta), path ("${prefix}.irmsd"), emit: irmsd
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${msa.baseName}"
    """
    export TEMP='./'

    t_coffee -other_pg irmsd \
        $msa \
        $args \
        -template_file $template > ${prefix}.irmsd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | awk '{gsub("Version_", ""); print \$3}')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${msa.baseName}"
    """
    touch ${prefix}.irmsd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | awk '{gsub("Version_", ""); print \$3}')
    END_VERSIONS
    """
}
