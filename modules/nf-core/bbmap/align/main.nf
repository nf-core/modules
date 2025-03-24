process BBMAP_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
    tuple val(meta), path(fastq)
    path ref

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    input = meta.single_end ? "in=${fastq}" : "in=${fastq[0]} in2=${fastq[1]}"

    // Set the db variable to reflect the three possible types of reference input: 1) directory
    // named 'ref', 2) directory named something else (containing a 'ref' subdir) or 3) a sequence
    // file in fasta format
    if ( ref.isDirectory() ) {
        if ( ref ==~ /(.\/)?ref\/?/ ) {
            db = ''
        } else {
            db = "path=${ref}"
        }
    } else {
        db = "ref=${ref}"
    }

    """
    bbmap.sh \\
        $db \\
        $input \\
        out=${prefix}.bam \\
        $args \\
        threads=$task.cpus \\
        -Xmx${task.memory.toGiga()}g \\
        &> ${prefix}.bbmap.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bbmap.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
