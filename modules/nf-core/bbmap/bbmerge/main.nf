process BBMAP_BBMERGE {
    tag "$meta.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.06--h92535d8_0':
        'biocontainers/bbmap:39.06--h92535d8_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_merged.fastq"), emit: merged
    tuple val(meta), path("*_unmerged.fastq"), emit: unmerged
    tuple val(meta), path("*_ihist.txt"), emit: ihist
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    bbmerge.sh \\
        -Xmx\$maxmem \\
        in1=${reads[0]} \\
        in2=${reads[1]} \\
        out=${prefix}_merged.fastq \\
        outu1=${prefix}_1_unmerged.fastq \\
        outu2=${prefix}_2_unmerged.fastq \\
        ihist=${prefix}_ihist.txt \\
        $args \\
        &> ${prefix}.bbmerge.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_merged.fastq
    touch ${prefix}_1_unmerged.fastq
    touch ${prefix}_2_unmerged.fastq
    touch ${prefix}_ihist.txt
    touch ${prefix}.bbmerge.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
