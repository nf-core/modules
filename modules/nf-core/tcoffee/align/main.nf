process TCOFFEE_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta) ,  path(fasta)
    tuple val(meta2),  path(tree)
    tuple val(meta3),  path(template), path(accessory_informations)

    output:
    tuple val (meta), path ("*.aln"), emit: msa
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """    
    t_coffee -seq ${fasta} \
        $args \
        -outfile ${prefix}.aln


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t_coffee: \$( t_coffee -version | sed 's/.*(Version_\\(.*\\)).*/\\1/' )
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t_coffee: \$( t_coffee -version | sed 's/.*(Version_\\(.*\\)).*/\\1/' )
    END_VERSIONS
    """
}
