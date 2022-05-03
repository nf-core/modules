process DRAGONFLYE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::dragonflye=1.0.11" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dragonflye:1.0.11--hdfd78af_0' :
        'quay.io/biocontainers/dragonflye:1.0.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("contigs.fa")                                        , emit: contigs
    tuple val(meta), path("dragonflye.log")                                    , emit: log
    tuple val(meta), path("{flye,miniasm,raven}.fasta")                        , emit: raw_contigs
    tuple val(meta), path("{miniasm,raven}-unpolished.gfa"), optional:true     , emit: gfa
    tuple val(meta), path("flye-info.txt"), optional:true                      , emit: txt
    path "versions.yml"                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def memory = task.memory.toGiga()
    """
    dragonflye \\
        --reads ${reads} \\
        $args \\
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
