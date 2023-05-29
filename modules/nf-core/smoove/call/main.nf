process SMOOVE_CALL {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::smoove=0.2.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/smoove:0.2.8--h9ee0642_1' :
        'biocontainers/smoove:0.2.8--h9ee0642_1' }"

    input:
    tuple val(meta), path(input), path(index), path(exclude_beds)
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

    def exclude = exclude_beds ? "--exclude ${exclude_beds}" : ""
    """
    smoove call \\
        ${args} \\
        --outdir . \\
        --name ${prefix} \\
        --fasta ${fasta} \\
        ${exclude} \\
        --processes ${task.cpus} \\
        ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smoove: \$(echo \$(smoove -v) | sed 's/^.*version: //; s/ .*\$//' )
    END_VERSIONS
    """
}
