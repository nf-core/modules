process TCOFFEE_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a76a981c07359a31ff55b9dc13bd3da5ce1909c1:84c8f17f1259b49e2f7783b95b7a89c6f2cb199e-0':
        'biocontainers/mulled-v2-a76a981c07359a31ff55b9dc13bd3da5ce1909c1:84c8f17f1259b49e2f7783b95b7a89c6f2cb199e-0' }"

    input:
    tuple val(meta) ,  path(fasta)
    tuple val(meta2),  path(tree)
    tuple val(meta3),  path(template), path(accessory_informations)
    val(compress)

    output:
    tuple val(meta), path("*.aln{.gz,}"), emit: alignment
    // in the args there might be the request to generate a lib file, so the following is an optional output
    tuple val(meta), path("*.*lib")      , emit: lib, optional : true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tree_args = tree ? "-usetree $tree" : ""
    def template_args = template ? "-template_file $template" : ""
    def write_output = compress ? " >(pigz -cp ${task.cpus} > ${prefix}.aln.gz)" : "> ${prefix}.aln"
    // using >() is necessary to preserve the tcoffee return value,
    // so nextflow knows to display an error when it failed
    """
    export TEMP='./'
    t_coffee -seq ${fasta} \
        $tree_args \
        $template_args \
        $args \
        -thread ${task.cpus} \
        -outfile stdout \
        $write_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | awk '{gsub("Version_", ""); print \$3}')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln${compress ? '.gz':''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | awk '{gsub("Version_", ""); print \$3}')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
}
