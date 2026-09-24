process NEXTDENOVO {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nextdenovo:2.5.2--py310h0ceaa1d_6' :
        'biocontainers/nextdenovo:2.5.2--py310h0ceaa1d_6' }"

    input:
    tuple val(meta), path(reads)
    path config

    output:
    tuple val(meta), path("*.fasta.gz"), emit: fasta
    tuple val(meta), path("*.stat")     , emit: stat
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "parallel_jobs = ${task.cpus}" >> conf.cfg
    cat $config >> conf.cfg
    echo ${reads} > input.fofn
    nextDenovo \\
        conf.cfg \\
        input.fofn \\

    gzip -c ./03.ctg_graph/nd.asm.fasta > ${prefix}.assembly.fasta.gz
  
    mv ./03.ctg_graph/nd.asm.fasta.stat ${prefix}.assembly_info.stat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$( nextDenovo --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo stub | gzip -c > ${prefix}.assembly.fasta.gz
    echo contig_1 > ${prefix}.assembly_info.stat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    \$( nextDenovo --version )
    END_VERSIONS
    """
}
