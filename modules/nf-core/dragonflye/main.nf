process DRAGONFLYE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dragonflye:1.1.2--hdfd78af_0' :
        'biocontainers/dragonflye:1.1.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(shortreads), path(longreads)

    output:
    tuple val(meta), path("*.fa")                                               , emit: contigs
    tuple val(meta), path("dragonflye.log")                                     , emit: log
    tuple val(meta), path("{flye,miniasm,raven}.fasta")                         , emit: raw_contigs
    tuple val(meta), path("{flye,miniasm,raven}-unpolished.gfa"), optional:true , emit: gfa
    tuple val(meta), path("flye-info.txt")                      , optional:true , emit: txt
    path "versions.yml"                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
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
}
