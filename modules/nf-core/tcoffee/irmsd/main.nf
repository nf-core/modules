process TCOFFEE_IRMSD {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/t-coffee_pigz:91ac7e26b23bb246':
        'community.wave.seqera.io/library/t-coffee_pigz:7d1373a24f76afe6' }"

    input:
    tuple  val(meta),  path (msa)
    tuple  val(meta2), path(template), path(structures)

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

    if [[ \$(basename $msa) == *.gz ]] ; then
        unpigz -f $msa
    fi

    t_coffee -other_pg irmsd \
        \$(basename $msa .gz) \
        $args \
        -template_file $template > ${prefix}.irmsd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | awk '{gsub("Version_", ""); print \$3}')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${msa.baseName}"
    """
    # Otherwise, tcoffee will crash when calling its version
    export TEMP='./'
    touch ${prefix}.irmsd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | awk '{gsub("Version_", ""); print \$3}')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
}
