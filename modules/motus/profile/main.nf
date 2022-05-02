process MOTUS_PROFILE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::motus=3.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/motus:3.0.1--pyhdfd78af_0':
        'quay.io/biocontainers/motus:3.0.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path db
    path bam

    output:
    tuple val(meta), path("*.out"), emit: out
    tuple val(meta), path("*.bam"), optional: true, emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def inputs = "$reads[0]".toLowerCase().endsWith('.bam') ?
                    "-i ${reads}" :
                    meta.single_end ?
                        "-s $reads" : "-f ${reads[0]} -r ${reads[1]}"
    def refdb = db ? "-db ${db}" : ""
    def intermediateBam = bam ? "-I $bam" : ""
    """
    motus profile \\
        $args \\
        $inputs \\
        $refdb \\
        $intermediateBam \\
        -t $task.cpus \\
        -n $prefix \\
        -o ${prefix}.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mOTUs: \$(echo \$(motus -h 2>&1) | sed 's/^.*Version: //; s/References.*\$//')
    END_VERSIONS
    """
}
