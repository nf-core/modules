process VGAN_HAPLOCART {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vgan:3.0.0--h9ee0642_0':
        'community.wave.seqera.io/library/vgan:3.1.0--cceaf026719e4db9' }"

    input:
    tuple val(meta), path(reads)
    val(is_interleaved)

    output:
    tuple val(meta), path("*.output.txt") , emit: txt
    tuple val(meta), path("*.posterior.txt")   , emit: posterior
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_args = (meta.single_end || is_interleaved ) ? "-fq1 ${reads}" : "-fq1 ${reads[0]} -fq2 ${reads[1]}"
    def interleaved = is_interleaved ? "-i" : ""
    """
    mkdir tmp
    vgan haplocart \\
        $args \\
        -t $task.cpus \\
        $reads_args \\
        $interleaved \\
        -o ${prefix}.output.txt \\
        -pf ${prefix}.posterior.txt \\
        -z tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vgan: \$(vgan version 2>&1 | sed -e "s/vgan version //g;s/ (Mela)//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.output.txt
    touch ${prefix}.posterior.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vgan: \$(vgan version 2>&1 | sed -e "s/vgan version //g;s/ (Mela)//g")
    END_VERSIONS
    """

}
