process MOTUS_PROFILE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::motus=3.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/motus:3.0.1--pyhdfd78af_0' :
        'quay.io/biocontainers/motus:3.0.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("${prefix}_motus.out")         , emit: motus_out
    tuple val(meta), path("${prefix}_motus.out.krona")   , emit: krona_input
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def inputs = meta.single_end ? "-s $reads" : "-f ${reads[0]} -r ${reads[1]}"
    """
    motus profile \\
        $inputs \\
        -t $task.cpus \\
        -n ${prefix} \\
        -db ${db} \\
        $args \\
        -I ${prefix}.bam \\
        -o ${prefix}_motus.out
    # for krona
    motus profile \\
        -i ${prefix}.bam \\
        -t $task.cpus \\
        -n ${prefix} \\
        -db ${db} \\
        -C parenthesis \\
        -o ${prefix}_motus.out.krona
    # clean up
    rm ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mOTUs: \$(grep motus ${db}/db_mOTU_versions | sed 's/motus\\t//g')
    END_VERSIONS
    """
}
