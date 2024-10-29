process NEXTGENMAP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nextgenmap%3A0.5.5--hc9558a2_4' :
        'biocontainers/nextgenmap:0.5.5--hc9558a2_4' }"

    input:
    tuple val(meta), path(reads)
    path(fasta)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus

    if(meta.single_end){
        """
        ngm \\
            -r $fasta \\
            -q $reads \\
            -t $threads \\
            --bam \\
            -o ${prefix}.bam \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            NextGenMap: \$(ngm 2>&1 | head -1 | grep -o -E '[[:digit:]]+.[[:digit:]]+.[[:digit:]]+')
        END_VERSIONS
        """
    } else{
        """
        ngm \\
            -r $fasta \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -t $threads \\
            --bam \\
            -o ${prefix}.bam \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            NextGenMap: \$(ngm 2>&1 | head -1 | grep -o -E '[[:digit:]]+.[[:digit:]]+.[[:digit:]]+')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NextGenMap: \$(ngm 2>&1 | head -1 | grep -o -E '[[:digit:]]+.[[:digit:]]+.[[:digit:]]+')
    END_VERSIONS
    """
}
