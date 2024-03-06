process BBMAP_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0':
        'biocontainers/bbmap:39.01--h5c4e2a8_0' }"

    input:
    tuple val(meta), path(fastq)
    
    output:
    tuple val(meta), path("*_merged.fastq"), emit: merged_fastq
    tuple val(meta), path("*_unmerged_R1.fastq"), emit: unmerged_fastq_1
    tuple val(meta), path("*_unmerged_R2.fastq"), emit: unmerged_fastq_2
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$meta.id"

    input = meta.single_end ? "in=${fastq}" : "in=${fastq[0]} in2=${fastq[1]}"

    """
    bbmerge.sh \\
        $input \\
	out=${prefix}_merged.fastq \\
	outu=${prefix}_unmerged_R1.fastq \\
	outu2=${prefix}_unmerged_R2.fastq \\
        $args \\
        threads=$task.cpus \\
        -Xmx${task.memory.toGiga()}g \\
	&> ${prefix}.bbmerge.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
