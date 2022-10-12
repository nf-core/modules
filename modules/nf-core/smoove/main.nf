process SMOOVE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::smoove=0.2.8" : null)
    container "brentp/smoove"

    input:
    tuple val(meta), path(input), path(index)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    smoove call \\
        ${args} \\
        --outdir . \\
        --name ${prefix} \\
        --fasta ${fasta} \\
        -p $task.cpus \\
        ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smoove: 0.2.8
    END_VERSIONS
    """
}
