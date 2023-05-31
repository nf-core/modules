process SPACERANGER_COUNT {
    tag "$meta.id"
    label 'process_high'

    // TODO push to nf-core docker
    container "ghcr.io/grst/spaceranger:2.0.1"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "SPACERANGER_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }


    input:
    tuple val(meta), path(fastq_dir), path(image), path(alignment)
    path(reference)
    path(probeset)

    output:
    tuple val(meta), path("**/outs/**", emit: outs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    spaceranger count \\
        --id="${meta.id}" \\
        --sample="${meta.id}" \\
        --fastqs="${fastq_dir}" \\
        --image="${image}" \\
        --slide="${meta.slide}" \\
        --area="${meta.area}" \\
        --transcriptome="${reference}" \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        $args
 
    samtools \\
        sort \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.bam \\
        -T $prefix \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
