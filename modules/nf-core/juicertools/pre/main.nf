process JUICERTOOLS_PRE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a2/a268a257cdea987bd60f7717686134f1a3c949e2ae268284642f1ce5a0434289/data' :
        'community.wave.seqera.io/library/juicertools_openjdk:fe58dd49794d6603' }"

    input:
    tuple val(meta) , path(pairs)
    tuple val(meta2), val(genome_id), path(chromsizes)

    output:
    tuple val(meta), path("*.hic"), emit: hic
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if(genome_id && chromsizes) {
        log.error("Error: both genome_id and chromsizes provided to juicertools/pre! Only one of these may be specified.")
    }
    if(!genome_id && !chromsizes) {
        log.error("Error: neither genome_id nor chromsizes provided to juicertools/pre!")
    }
    def args     = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    input_genome = genome_id       ?: chromsizes
    """
    export JAVA_OPTS='"-Xms${task.memory.toMega()/4}m" "-Xmx${task.memory.toGiga()}g"'

    juicer_tools pre \\
        --threads ${task.cpus} \\
        ${args} \\
        ${pairs} \\
        ${prefix}.hic \\
        ${input_genome}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        juicer_tools: \$(juicer_tools -V | grep "Juicer Tools Version" | sed 's/Juicer Tools Version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hic

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        juicer_tools: \$(juicer_tools -V | grep "Juicer Tools Version" | sed 's/Juicer Tools Version //')
    END_VERSIONS
    """
}
