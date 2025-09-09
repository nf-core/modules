process BISMARK_METHYLATIONEXTRACTOR {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe2a9d58209a38df5c99615a41d9ed6e8e546380d04c176e076e107010819a72/data' :
        'community.wave.seqera.io/library/bismark:0.25.0--95ba99b483e2eaf9' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bedGraph.gz")         , emit: bedgraph
    tuple val(meta), path("*.txt.gz")              , emit: methylation_calls
    tuple val(meta), path("*.cov.gz")              , emit: coverage
    tuple val(meta), path("*_splitting_report.txt"), emit: report
    tuple val(meta), path("*.M-bias.txt")          , emit: mbias
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Assign sensible numbers for multicore and buffer_size based on bismark docs
    if(!args.contains('--multicore') && task.cpus >= 6){
        args += " --multicore ${(task.cpus / 3) as int}"
    }
    // Only set buffer_size when there are more than 6.GB of memory available
    if(!args.contains('--buffer_size') && task.memory?.giga > 6){
        args += " --buffer_size ${task.memory.giga - 2}G"
    }

    def seqtype  = meta.single_end ? '-s' : '-p'
    """
    bismark_methylation_extractor \\
        ${bam} \\
        --bedGraph \\
        --counts \\
        --gzip \\
        --report \\
        ${seqtype} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bedGraph.gz
    touch ${prefix}.txt.gz
    touch ${prefix}.cov.gz
    touch ${prefix}_splitting_report.txt
    touch ${prefix}.M-bias.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
