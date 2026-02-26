process DRAGONFLYE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dragonflye:1.2.1--hdfd78af_0' :
        'biocontainers/dragonflye:1.2.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(shortreads), path(longreads)

    output:
    tuple val(meta), path("*.fa")                               , emit: contigs
    tuple val(meta), path("dragonflye.log")                     , emit: log
    tuple val(meta), path("{flye,miniasm,raven}.fasta")         , emit: raw_contigs
    tuple val(meta), path("{flye,miniasm,raven}-unpolished.gfa"), emit: gfa, optional:true
    tuple val(meta), path("flye-info.txt")                      , emit: txt, optional:true
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def memory  = task.memory.toGiga()
    def shortreads_polishing = shortreads ? "--R1 ${shortreads[0]} --R2 ${shortreads[1]}" : ''
    """
    dragonflye \\
        --reads ${longreads} \\
        $shortreads_polishing \\
        $args \\
        --prefix ${prefix} \\
        --cpus $task.cpus \\
        --ram $memory \\
        --outdir ./ \\
        --force

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragonflye: \$(dragonflye --version 2>&1 | sed 's/^.*dragonflye //' )
    END_VERSIONS
    """

    stub:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def args      = task.ext.args   ?: ''
    def args_list = args.tokenize()

    assembler_name = args_list.contains("--assembler") ? args_list[args_list.indexOf('--assembler')+1] : 'flye'

    """
    touch ${prefix}.fa
    touch ${prefix}.reoriented.fa
    touch dragonflye.log
    touch ${assembler_name}.fasta
    touch ${assembler_name}-unpolished.gfa
    if [ "${assembler_name}" == "flye" ]; then
        touch flye-info.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragonflye: \$(dragonflye --version 2>&1 | sed 's/^.*dragonflye //' )
    END_VERSIONS
    """
}
