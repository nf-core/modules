process LEARNMSA_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-741e0da5cf2d6d964f559672e2908c2111cbb46b:4930edd009376542543bfd2e20008bb1ae58f841-0' :
        'biocontainers/mulled-v2-741e0da5cf2d6d964f559672e2908c2111cbb46b:4930edd009376542543bfd2e20008bb1ae58f841-0' }"

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
    learnMSA \\
        $args \\
        -i <(unpigz -cdf $fasta) \\
        -o $write_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        learnmsa: \$(learnMSA -h | grep 'version' | awk -F 'version ' '{print \$2}' | awk '{print \$1}' | sed 's/)//g')
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
        learnmsa: \$(learnMSA -h | grep 'version' | awk -F 'version ' '{print \$2}' | awk '{print \$1}' | sed 's/)//g')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
}
