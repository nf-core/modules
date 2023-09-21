process TCOFFEE_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::t-coffee=13.45.0.4846264"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/t-coffee:13.45.0.4846264--hc57179f_5':
        'biocontainers/t-coffee:13.45.0.4846264--hc57179f_5' }"

    input:
    tuple val(meta) ,  path(fasta)
    tuple val(meta2),  path(tree)
    tuple val(meta3),  path(template), path(accessory_informations)

    output:
    tuple val (meta), path ("*.aln"), emit: msa
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """    
    t_coffee -seq ${fasta} \
        $args \
        -outfile ${prefix}.aln


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | sed 's/.*(Version_\\(.*\\)).*/\\1/' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | sed 's/.*(Version_\\(.*\\)).*/\\1/' )
    END_VERSIONS
    """
}
