process KALIGN_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5cd0277547c6b33133225c8ce14c0cf2a4396ea2:0a70b6d89a3e06fbdc4a735461e8b98ff32ee5de-0':
        'biocontainers/mulled-v2-5cd0277547c6b33133225c8ce14c0cf2a4396ea2:0a70b6d89a3e06fbdc4a735461e8b98ff32ee5de-0' }"

    input:
    tuple val(meta), path(fasta)
    val(compress)

    output:
    tuple val(meta), path("*.aln{.gz,}"), emit: alignment
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def write_output = compress ? ">(pigz -cp ${task.cpus} > ${prefix}.aln.gz)" : "${prefix}.aln"
    """
    error_handler() {
    exit_code=\$?
        if [ \$exit_code -eq 132 ]; then
            echo "KALIGN failed because is incompatible with some CPU types, see https://github.com/TimoLassmann/kalign/issues/46."
        else
            trap - ERR
            return \$exit_code
        fi
    }

    trap 'error_handler' ERR

    unpigz -cdf $fasta | \\
    kalign \\
        $args \\
        -n ${task.cpus} \\
        -o ${write_output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kalign: \$(echo \$(kalign -v) | sed 's/kalign //g' )
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln${compress ? '.gz' : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kalign : \$(echo \$(kalign -v) | sed 's/kalign //g' )
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
}
