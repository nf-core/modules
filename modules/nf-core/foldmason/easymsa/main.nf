process FOLDMASON_EASYMSA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/foldmason_pigz:97b3311addb0f4a7"

    input:
    tuple val(meta), path(pdbs)
    val(compress)

    output:
    tuple val(meta), path("${prefix}_3di.fa${compress ? '.gz' : ''}"), emit: msa_3di
    tuple val(meta), path("${prefix}_aa.fa${compress ? '.gz' : ''}") , emit: msa_aa
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    foldmason easy-msa \\
        $args \\
        --threads $task.cpus \\
        ${pdbs} \\
        ${prefix} \\
        tmp

    if ${compress}; then
        pigz -p ${task.cpus} ${prefix}_3di.fa
        pigz -p ${task.cpus} ${prefix}_aa.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldmason: \$(foldmason | grep "foldmason Version:" | cut -d":" -f 2 | awk '{\$1=\$1;print}')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ""  ${compress ? '| gzip' : ''} > ${prefix}_3di.fa${compress ? '.gz' : ''}
    echo ""  ${compress ? '| gzip' : ''} > ${prefix}_aa.fa${compress ? '.gz' : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldmason: \$(foldmason | grep "foldmason Version:" | cut -d":" -f 2 | awk '{\$1=\$1;print}')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
}
