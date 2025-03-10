process BBMAP_BBNORM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bbmap_samtools_pigz:2a066f0214cc5eb0' :
        'community.wave.seqera.io/library/bbmap_samtools_pigz:79703e935236b43b' }"
    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    tuple val(meta), path("*.log")     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    input  = meta.single_end ? "in=${fastq.join(',')}" : "in=${fastq[0]} in2=${fastq[1]}"
    output = meta.single_end ? "out=${prefix}.fastq.gz" : "out1=${prefix}_1.nm.fastq.gz out2=${prefix}_2.nm.fastq.gz"

    """
    bbnorm.sh \\
        $input \\
        $output \\
        $args \\
        threads=$task.cpus \\
        -Xmx${Math.round(Math.max(1, Math.floor(task.memory.toGiga() * 0.95)))}g \\
        &> ${prefix}.bbnorm.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
